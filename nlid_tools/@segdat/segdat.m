classdef segdat<nldat
    %Segmented nldat object
    properties
        segInfo ={};
        parameterSet

    end
    methods
        function S = segdat (a,varargin)
            %segdat specific parameters
            S.parameterSet=param('paramName','onsetPointer','paramDefault',1, ...
                'paramHelp','a vector of segments onsets (sample index)', ...
                'paramType','number');
            S.parameterSet(2)=param('paramName','segLength','paramDefault',0, ...
                'paramHelp','a vector of segments lengths (number of samples in each segment)', ...
                'paramType','number');
            if nargin==0;
                return
            elseif isa(a,'double')
                S=segdat.nl2seg(nldat(a),inputname(1));
                if nargin > 1,
                    set (S,varargin{:});
                end
            elseif isa(a,'segdat')
                S = nlmkobj(a,varargin{:});
            elseif isa(a,'nldat')
                S=segdat.nl2seg(a,inputname(1));
            else
                set (S,{a varargin{:}});
            end

        end

        function Y = resampleSeg(X,fs,chan)
            % Definition: resamples the data set of a segdat object, one segment at a
            % time. Resampling is based on a specified channel containing the time
            % values, or the domain properties of the object itself. This function
            % resamples the entire data set.
            % Inputs:
            %   X = signal to be resampled (as a segdat object)
            %   fs = desired sampling frequency (Hz)
            %   chan = channel number containing the time values (optional -- if
            %   not specified, a vector of times will be created using the domain
            %   starts of each segment, the domain increments and the segment lengths)
            % Output:
            %   Y = resampled signal

            % Extract pertinent information
            nSeg = segCount(X);
            domainIncr = X.domainIncr;
            comment = X.comment;
            chanNames = X.chanNames;
            chanUnits = X.chanUnits;

            % Read input arguments
            ts = 1/fs; % Sample time (sec.)

            for i = 1:nSeg

                seg_i = segGet(X,i); % Extract segment
                seg_data = seg_i.dataSet; % Copy the data into a manipulable variable

                % Establish the time vector
                if ~exist('chan','var')
                    % If no channel number specified, create a vector of uniformly
                    % sampled times based on the start time and length of the segment
                    time = domain(seg_i); % Create an array to store the time values of the segment
                else
                    time = seg_data(:,chan); % If a channel number was specified, use the times in that channel for resampling
                end

                % Check the uniformity of the segment data
                dt = diff(time);

                if range(dt) > 1e-5
                    display(['Note: segment ',num2str(i),' is non-uniformly sampled. Correcting data...']);
                else
                    display(['Note: segment ',num2str(i),' is uniformly sampled.']);
                end

                % Pad the front of the segment
                nPad = 10; % Number of samples to be padded in front and back of segment
                tFront = [time(1)-nPad*domainIncr:domainIncr:time(1)-domainIncr]';
                padFront = zeros(nPad,size(seg_data,2));
                for k = 1:nPad
                    padFront(k,:) = seg_data(1,:);
                end

                % Pad the end of the segment
                tBack = [time(end)+domainIncr:domainIncr:time(end)+nPad*domainIncr]';
                padBack = zeros(nPad,size(seg_data,2));
                for k = 1:nPad
                    padBack(k,:) = seg_data(end,:);
                end

                % Apply the padding to both the time vector and the segment data set
                time_ext = [tFront;time;tBack];
                seg_data_ext = [padFront;seg_data;padBack];

                % If necessary, correct the time stamps in the segement data set
                if exist('chan','var')
                    seg_data_ext(:,chan) = time_ext;
                end

                % Resample the entire data set of the segment plus the padding
                [new_data,new_time] = resample(seg_data_ext,time_ext,fs);

                % Smooth the resampled signal by applying a moving average filter
                nMA = 2; % Window size of moving average filter (in sample points)
                B = 1/nMA*ones(nMA,1); % Filter coefficients
                new_data = filter(B,1,new_data);
                % new_data = medfilt1(new_data);

                % Round the time stamps to integer multiples of the new sampling time
                new_time = ts*round(new_time/ts);

                % If the time vector channel was specified, replace the channel values
                % with the new time vector
                if exist('chan','var')
                    new_data(:,chan) = new_time;
                end

                % Find the first sample at/after the original start time of the segment
                % post-resampling

                trimF=min(find((new_time>=time(1))))
                % Trim the beginning of the data set
                %     new_data(1:(trimF-1),:) = [];
                %     new_time(1:(trimF-1),:) = [];

                % Find the first sample at/before the original end time of the segment
                % post-resampling

                trimB=min(find(new_time>=max(time)))

                % Extract the portion of segment between [trimF,trimB]
                new_data=new_data(trimF:trimB,:);
                new_time=new_time(trimF:trimB,:);

                % Overwrite data in extracted segment
                seg_i.dataSet = new_data;

                % Correct segment properties post-resampling
                set(seg_i,'domainStart',new_time(1),'domainIncr',ts,'comment',comment,'chanNames',chanNames,'chanUnits',chanUnits);

                % Create a new segdat object or choose to link the segment to
                % the current one
                if i == 1
                    Y = segdat(seg_i); % If this is the first segment, create the segdat object Y
                else
                    Y = cat(1,Y,seg_i); % Otherwise, link the segment to Y
                end

                display(['Segment ', num2str(i),' successfully resampled and saved.']);

            end

        end

        % Last updated: 27-05-2022

        function N=seg4domain (S, domainVal)
            % Return segment mumber associated with a domain Value
            dStart=S.domainStart;
            dEnd=domainEnd(S);
            delta=-S.domainIncr/2;
            N=find((domainVal-dStart)>=delta & (domainVal-dEnd)<=delta);
        end
        function C=cor(S,varargin)
            C=corirf(S,'cor',varargin{:});
        end
        function I=irf(S,varargin)
            I=corirf(S,'irf',varargin{:});
        end


        function e= segdat2eseq (S, domainStart)
            % e= segdat2eseq (S, domainStart)
            % domainstart for sequence data (seconds);
            e=eseq;
            nSeg=segCount(S);
            for iSeg=1:nSeg
                curSeg=segGet(S,iSeg);
                curDomain=domain(curSeg);
                curIdx=idx4domain(domainStart, S.domainIncr, curDomain);
                e(iSeg).domainStart=domainStart;
                e(iSeg).domainIncr=S.domainIncr;
                e(iSeg).startIdx=min(curIdx);
                e(iSeg).endIdx=max(curIdx);
                e(iSeg).nSamp=length(curSeg);
                e(iSeg).type='SEGDAT';
            end
        end



        function C = corirf(S, fncType,varargin)
            % segdatCor - overlaid correlation function for segdat objecs
            %   Detailed explanation goes here
            iLag=find(strcmp('nLags',varargin));
            if isempty(iLag)
                error('nlags must be specified for cor on segdat objects');
            else
                nLags=varargin{iLag+1};
            end
            nSeg=segCount(S);
            C={};
            iCnt=0;
            for iSeg=1:nSeg
                curSeg=segGet(S,iSeg);
                curLen=length(curSeg);
                if curLen>2*nLags+1
                    switch fncType
                        case 'cor'
                            cTemp=cor(curSeg, varargin);
                        case 'irf'
                            cTemp=irf(curSeg, varargin);
                        otherwise
                            error(['Bad fncType':  fncType])
                    end

                    iCnt=iCnt+1;
                    if iCnt==1,
                        C=cTemp;
                    else
                        C=C+cTemp;
                    end
                else
                    disp(['Segment ' num2str(iSeg) ' to0 short not included']);
                end
            end
            if iCnt>0
                C=C./iCnt;
                set(C,'comment',fncType);
            else
            end
        end

        function sDelay = delay (S, nDelay)
            % delay all segemnts of S by nDelay samples
            sDelay=S;
            nSeg=segCount(S);
            for iSeg=1:nSeg,
                curSeg=segGet(S,iSeg);
                curSegD=delay(curSeg,nDelay);
                curData=double(curSegD);
                if iSeg==1,
                    dataSet=curData;
                else
                    dataSet=cat(1,dataSet,curData);

                end
            end
            set(sDelay,'dataSet', dataSet);
        end

        function plot(S)
            [nSamp,nChan,nReal]=size(S);
            colors=colororder;
            numSegment = segCount(S);
            if (nChan==1 && nReal==1)

                % plot(S_nldat)
                hold on
                for i = 1 : numSegment
                    sSeg=segGet(S,i);
                    plot(sSeg)
                end
                hold off

            elseif nChan>1
                for iChan=1:nChan
                    subplot (nChan, 1,iChan)
                    hold on
                    for iSeg = 1 : numSegment
                        sSeg=segGet(S,iSeg);
                        plot(sSeg(:,iChan))
                        title(S.comment)
                    end
                    hold off
                end
            end
        end


        function sOut = segBreak (sIn, segNum, iBreak)
            %  sOut = segBreak (sIn, segNum, iBreak)
            %           break a segment into two pieces
            curSeg=segGet(sIn,segNum);
            dSeg=domain(curSeg);
            segLen=length(curSeg);
            newSegNum=segNum+1;
            oldDomainStart=sIn.domainStart;
            newDomainStart=vinsert(oldDomainStart,segNum+1,dSeg(iBreak));
            oldOnset=get(sIn,'onsetPointer');
            newOnsetPointer=vinsert(oldOnset, newSegNum,oldOnset(segNum)+iBreak-1);
            oldSegLen=get(sIn,'segLength');
            oldSegLen(segNum)=iBreak-1;
            newSegLength=vinsert(oldSegLen,newSegNum, segLen-iBreak+1);
            sOut=sIn;
            sOut.domainStart=newDomainStart;
            set(sOut,'onsetPointer',newOnsetPointer);
            set(sOut,'segLength',newSegLength);
        end

        function Z = segCat (Z1, Z2, varargin)
            % Concatonate segdat objects - only for dim 1.
            % In case of overlap the output is given by Z2
            if nargin>2,
                Z=segCat(Z1,Z2);
                for iArg=3:nargin,
                    jArg=iArg-2;
                    Z=segCat(Z,varargin{jArg});
                end
            else
                domainIncr=Z1.domainIncr;
                if Z2.domainIncr ~= domainIncr
                    error('DomainIncrements are not the same');
                end
                name1=inputname(1);
                name2=inputname(2);
                n1=nldat(Z1);
                d1=domain(n1);
                n2=nldat(Z2);
                d2=domain(n2);
                d=cat(1,d1,d2);
                dMin=min(d);
                dMax=max(d);
                iMax=idx4domain(dMin, domainIncr, dMax);
                zd=nan(iMax,1);
                ptr1=idx4domain(dMin,domainIncr,d1);
                zd(ptr1,:)=n1(:,:);
                ptr2=idx4domain(dMin,domainIncr,d2);
                iOverlap=intersect(ptr1,ptr2);
                comment=['segcat(' name1 ',' name2 ')'];
                if ~isempty (iOverlap),
                    overlapStart=min(iOverlap);
                    overlapEnd=max(iOverlap);
                    comment=['segcat:' name2 ' overlaps ' name1 ' from ' num2str(overlapStart) ...
                        ':' num2str(overlapEnd)];
                end
                zd(ptr2,:)=n2(:,:);
                set (n1,'domainStart',dMin,'dataSet',zd, 'comment' ,comment);
                Z=segdat(n1);
            end
        end

        function dEnd =domainEnd (sIn)
            % Returns end domain values for segments in a segdat object
            domainStart=sIn.domainStart;
            segLen=get(sIn,'segLength');
            dEnd=domainStart+(segLen-1)*sIn.domainIncr;
        end



        function Z = intersect (Z1, Z2)
            % Determine intersection between two segdat objects
            % Returns an empty variable if the intersection is emmpty
            domainIncr=Z1.domainIncr;
            if Z2.domainIncr ~= domainIncr
                error('DomainIncrements are not the same');
            end
            name1=inputname(1);
            name2=inputname(2);
            n1=nldat(Z1);
            iVALID1=find(~isnan(double(n1(:,1))));
            d1=domain(n1);
            d1=d1(iVALID1);
            n2=nldat(Z2);
            iVALID2=find(~isnan(double(n2(:,1))));
            d2=domain(n2);
            d2=d2(iVALID2);
            iIntersect=ismember(d1,d2);
            if ~any(iIntersect),
                Z=[];
                return;
            end
            dCommon=d1(iIntersect);
            ptr1=idx4domain(n1.domainStart,domainIncr,dCommon);
            ptr2=idx4domain(n2.domainStart,domainIncr,dCommon);
            n1New=n1(ptr1,:);
            n2New=n2(ptr2,:);
            Z=cat(2,n1New,n2New);
        end



        function n = segCount(S)
            % reuturns number of segments in a segdat object
            onsetPointer=get(S,'onsetPointer');
            n=length(onsetPointer);
        end

        function N = segGet(S,segNum)
            % Return segment of segdat object as a nldat object
            N=nldat;
            nSeg=segCount(S);
            if (segNum>nSeg)
                error('too many segments')
            end
            dataSet=double(S);
            segLen=get(S,'segLength');
            onsetPointer=get(S,'onsetPointer');
            sStart=onsetPointer(segNum);
            sEnd=onsetPointer(segNum)+segLen(segNum)-1;
            N.dataSet=dataSet(sStart:sEnd,:);
            set(N, 'chanNames', S.chanNames,'chanUnits',S.chanNames, 'domainIncr',S.domainIncr, ...
                'domainName',S.domainName);
            N.comment=(['Segment ' num2str(segNum) '  of ' S.comment ]);
            N.domainStart=S.domainStart(segNum);
        end

        function segPlot (S)
            % plot segdat object
            nSeg= segCount(S);
            [nSamp,nChan]=size(S);
            for iChan=1:nChan,
                for iSeg=1:nSeg,
                    sTemp=segGet(S,iSeg);
                    line(domain(sTemp),sTemp(:,iChan));
                end
                xlabel(S.domainName);
                ylabel(S.chanNames{iChan});
            end
        end


        function d =domain(S)
            % Overlaid domain function for segdat objects;
            nSeg=segCount(S);
            for iSeg=1:nSeg;
                curSeg=segGet(S,iSeg);
                if iSeg==1;
                    d=domain(curSeg);
                else
                    d=cat(1,d,domain(curSeg));
                end
            end
        end





        function out = decimate(data,decimation_ratio)
            errorcheck(data);
            [~,nchan,~] =size(data);
            onsetpointer = get(data,'onsetPointer');
            seglength = get(data,'segLength');
            %             if size(seglength,2)==1
            %             elseif  any(~(mean(seglength')==seglength(:,1)'))
            %                 error('All data channels must have equal segment lengths')
            %             end
            %onsetpointer_new = zeros(size(onsetpointer));
            %seglength_new = zeros(size(onsetpointer));
            endpointer = onsetpointer + seglength - 1;
            ts = get(data,'domainIncr');
            out = data;
            set(out,'domainIncr',ts*decimation_ratio);
            dataset = get(data,'dataSet');
            %d = zeros(ceil(sum(seglength(:,1))/decimation_ratio),nchan);
            for i = 1: nchan
                pointer = 1;
                jOut=0;
                for j = 1 : length(onsetpointer)
                    curSeg=dataset(onsetpointer(j):endpointer(j),i);
                    if length(curSeg)<24,
                        disp(['Segment ' num2str(j)  ' too short. Dropping']);
                    else
                        jOut=jOut+1;
                        d_temp = decimate(curSeg,decimation_ratio);
                        onsetpointer_new(jOut) = pointer;
                        seglength_new(jOut) = length(d_temp);
                        d(onsetpointer_new(jOut):onsetpointer_new(jOut)+seglength_new(jOut)-1,i) = d_temp;
                        pointer = pointer + seglength_new(jOut);
                    end
                end
            end
            set(out,'dataSet',d,'onsetPointer',onsetpointer_new,'segLength',seglength_new);
        end

        function out = ddt(data)
            errorcheck(data);
            out=data;
            newDataSet=[];
            [nSamp,nChan,nReal]=size(data);
            for iChan=1:nChan,
                curData=data;
                curData.chanNames=data.chanNames{iChan};
                curData.dataSet=data.dataSet(:,iChan);
                nSeg=segCount(curData);
                xdt=nldat;
                for iSeg=1:nSeg
                    xTmp=segGet(curData,iSeg);
                    xdt=ddt(xTmp);
                    if iSeg==1,
                        curOut=segdat(xdt);
                    else
                        curOut=segCat(curOut,xdt);
                    end
                end
                newDataSet=cat(2, newDataSet,curOut.dataSet);
            end
            out.dataSet=newDataSet;
            out.comment='ddt';
        end

        function sOut =filter (sIn, varargin)
            % length(epochIn)+1) - Filter a segdat object
            % Filters each segmenet separately.
            % sIn - segdat object
            % vrargin -
            %     I - impulse repsonse of filter
            %     B,A - filter characteristics
            arg1=varargin{1};
            nSeg=segCount(sIn);
            for iSeg=1:nSeg,
                curSeg=segGet(sIn,iSeg);
                if isa(arg1,'irf')
                    curFilt=nlsim(arg1,curSeg);
                else
                    curFilt=filter(curSeg,varargin);
                end
                if iSeg==1,
                    sOut=segdat(curFilt);
                else
                    sOut=segCat(sOut,curFilt);
                end
            end

        end

        function [ fullEpochs, shortEpochs]= getEpochs ( S, epochLen)
            % [ fullEpochs, partialEpochs]= getEpochs ( S, epochLength fullEpochs=segdat; set(fullEpochs,'domainIncr',S.domainIncr);)
            % Break a segdat object into contiguous epochs of length =epochLen
            % fullEpochs - segdat object containing epochs of length= epochLen
            % shortEpochs - segdat object containing epochs of length< epochLen
            % length=epochLen
            fullEpochs=[];
            shortEpochs=[];


            nSeg=segCount(S);
            for iSeg=1:nSeg
                curSeg=segGet(S,iSeg);
                curSegLen=length(curSeg);
                for iStart=1:epochLen:curSegLen
                    iEnd= iStart+epochLen-1;
                    if iEnd>curSegLen
                        shortEpochs=segdat.helperCat(shortEpochs, curSeg(iStart:curSegLen,:));
                    else
                        fullEpochs=segdat.helperCat(fullEpochs, curSeg(iStart:iEnd,:));
                    end
                end

            end
        end

        function sOut = medfilt1 (sIn, varargin)
            nSeg=segCount(sIn);
            for iSeg=1:nSeg,
                curSeg=segGet(sIn,iSeg);
                curFilt=medfilt1(curSeg,varargin);
                if iSeg==1,
                    sOut=segdat(curFilt);
                else
                    sOut=segCat(sOut,curFilt);
                end
            end
        end


        function out = nldat(data)
            % returns concatonated segments as a nldat object with separations filled with
            % nan for missing data

            errorcheck(data);
            [~,nchan,~] =size(data);
            onsetPointer = get(data,'onsetPointer');
            seglength = get(data,'segLength');
            nSeg=segCount(data);
            [nSamp,nChan,nReal]=size(data);
            domainIncr=data.domainIncr;
            newDomain=nan;
            d=domain(data);
            domainStart=min(d);
            domainEnd=max(d);
            idx=idx4domain(domainStart,domainIncr, domainEnd);
            outData=nan(idx,nChan);
            for iSeg=1:nSeg,
                curSeg=segGet(data,iSeg);
                curIdx=idx4domain(domainStart,domainIncr, domain(curSeg));
                curData=double(curSeg);
                outData(curIdx,:)=curData(:,:);
            end
            out = nldat(outData,'chanNames',data.chanNames,'chanUnits',data.chanUnits, ...
                'domainIncr',data.domainIncr,'domainName',data.domainName,'domainStart',domainStart, ...
                'domainValues',nan,'comment',data.comment);
        end


        function v = vaf(data1,data2)
            [~,nchan1,~] =size(data1);
            [~,nchan2,~] =size(data2);
            if (nchan1>1) || (nchan2>1)
                error('VAF not supported for multiple channel data...')
            end
            d1 = nldat(data1.dataSet);
            d2 = nldat(data2.dataSet);
            v = vaf(d1,d2);
        end


        function errorcheck(data)
            dataset = get(data,'dataSet');
            [~,nchan,~] = size(data);
            onsetpointer = get(data,'onsetPointer');
            seglength = get(data,'segLength');
            if ~(size(onsetpointer) == size(seglength))
                error('onsetPointer and segLength must have equal dimension')
            end
            if ~isempty(find(seglength<0, 1))
                error('At least one segment does not have a correct length')
            end
            endpointer = onsetpointer + seglength - 1;
            if ~isempty(find(endpointer>size(dataset,1), 1))
                error('At least segment exceeds the data range')
            end
            %% Suppress checking for overlap since this is now OK
            %             for k = 1: nchan
            %                 for i = 1 : length(onsetpointer)
            %                     for j = 1 : length(onsetpointer)
            %                         if i == j
            %                             break;
            %                         elseif (onsetpointer(j)<=onsetpointer(i))&&(onsetpointer(i)<=endpointer(j))
            %                             % warning(['In channel ',num2str(k),', segment ',num2str(i),' & segment',num2str(j),' overlap.'])
            %                         end
            %                     end
            %                 end
            %
            %             end

        end


        function sCat = cat (DIM,S1,S2)
            % overlaid segdat function for segdat objects
            if isa(S2,'nldat')
                S2=segdat(S2);
            elseif ~isa(S2,'segdat')
                error('objects to be concatonated must both be segdat')
            end

            if DIM==1,
                % Concatonate two segdat objects
                [nSamp,nChan,nReal]=size(S1);
                domainIncr=S1.domainIncr;
                d1=domain(S1);
                s1Data=double(S1);
                d2=domain(S2);
                s2Data=double(S2);
                catVector =categorical;
                d2=domain(S2);
                domainStart=min(min(d1),min(d2));
                domainEnd=max(max(d1),max(d2));
                % Handle overlapping segments
                idxMax=idx4domain(domainStart, domainIncr, domainEnd);
                domainVector=domainStart:domainIncr:domainEnd;
                valueVector=zeros(length(domainVector),nChan);
                catVector(1:idxMax)='0';  % Initialize to nothing
                idx1=idx4domain(domainStart, domainIncr, d1); % s1 Values
                valueVector(idx1,:)=s1Data;
                catVector(idx1)='1';
                idx2=idx4domain(domainStart, domainIncr, d2);  % S2 values
                valueVector(idx2,:)=s2Data;
                catVector(idx2)='2';
                idxIntersect=intersect(idx1,idx2);
                if ~isempty(idxIntersect)  % Interesction
                    catVector(idxIntersect)='3';
                end
                e=eseq(catVector, domainStart,domainIncr); % Event sequences indicating what to output.
                eDomain=domain(e);
                % Generate concatonated segdat
                segNum=0;
                domainStart=[];
                onsetPointer=[];
                segInfo={};
                segLen=[];
                dataSet=[];
                for i=1:length(e),
                    if e(i).type ~='0'
                        segNum=segNum+1;
                        domainStart(segNum)=domainVector(e(i).startIdx);
                        iStart=e(i).startIdx;
                        iEnd=e(i).endIdx;
                        dataSet=cat(1, dataSet, valueVector(iStart:iEnd,:));
                        segLen= [ segLen e(i).nSamp];
                        eventStart=min(eDomain{i});
                        if segNum==1
                            onsetPointer(1)=1;
                        else
                            onsetPointer(segNum)=onsetPointer(segNum-1)+segLen(segNum-1);
                        end
                        switch e(i).type
                            case '1'
                                curSegNum=seg4domain(S1,eventStart);
                                segInfo{segNum}=S1.segInfo{curSegNum};
                            case '2'
                                curSegNum=seg4domain(S2,eventStart);
                                segInfo{segNum}=S2.segInfo{curSegNum};
                            case '3'
                                curSegNum=seg4domain(S2,eventStart);
                                segInfo{segNum}=S2.segInfo{curSegNum};
                        end
                    end
                end

                sCat=S1;
                sCat.domainStart=domainStart;
                sCat.dataSet=dataSet;
                set(sCat,'onsetPointer',onsetPointer,'segInfo',segInfo);
                set(sCat,'segLength',segLen);
            elseif DIM==2
                sCat=S1;
                sCat.dataSet=cat(2,S1.dataSet,S2.dataSet);
                sCat.chanNames= cat(2,S1.chanNames, S2.chanNames);
            else
                error ('Dimensions >2 not support');
            end
        end

        function out = subsref (N, S)

            nTemp=N;
            for i=1:length(S),
                if strcmp(S(i).type,'.'),
                    nTemp=get(nTemp,S(i).subs);
                elseif strcmp(S(i).type,'()')

                    oDom=domain(nTemp);
                    d=nTemp.dataSet;
                    dTemp=builtin('subsref', d, S(i));
                    nTemp.dataSet=dTemp;
                    sIndex = S(i).subs;
                    % Fix DomainStart and DomainValues if necessary
                    if isnumeric(sIndex{1}),
                        nDom = oDom(sIndex{1});
                        diffIndex=diff(sIndex{1});
                        nTemp.domainStart=min(nDom);
                        if ~isnan(nTemp.domainValues) | length(unique(diffIndex))>1,
                            nTemp.domainValues=nDom;
                        elseif max(diffIndex)>1,
                            nTemp.domainIncr=nTemp.domainIncr*max(diffIndex);
                        end
                    end
                    % Fix Channel Names
                    if length(sIndex)>1,
                        nTemp.chanNames = nTemp.chanNames(S(i).subs{2});
                    end
                    % Special handling for  segdat
                    if strcmp(S(i).subs(1),':')
                        return
                    end

                    onsetPointer=get(N,'onsetPointer');
                    segLength=get(N,'segLength')
                    nSeg=length(onsetPointer)
                    segnum4sample=[];
                    for iSeg=1:nSeg
                        segnum4sample=cat(1,segnum4sample,iSeg*ones(segLength(iSeg),1))
                    end
                    oldDomainValues=domain(N);
                    S(i).subs{2}=1;
                    segNumNew=builtin('subsref', segnum4sample, S(i));
                    newDomainValues=builtin('subsref', oldDomainValues, S(i));

                    e=eseq.cseq2eseq(categorical(segNumNew));
                    newOnsetPointer=[e.startIdx];
                    newSegLength=[e.nSamp];
                    newDomainStart(1,:)=newDomainValues(newOnsetPointer)

                    set(nTemp,'onsetPointer',newOnsetPointer);
                    set(nTemp,'segLength',newSegLength);
                    set(nTemp,'domainStart',newDomainStart);


                end
            end

            out = nTemp;
        end
    end




    methods (Static)

        function  epochOut = helperCat (epochIn, seg)
            if isempty (epochIn)
                epochOut=segdat(seg);
            else
                epochOut=epochIn;
                epochOut.domainStart=cat(2,epochIn.domainStart, seg.domainStart);
                epochOut.dataSet=cat(1,epochIn.dataSet, seg.dataSet);
                onsetPointer=get(epochIn,'onsetPointer');
                set(epochOut,'onsetPointer',cat(2,onsetPointer,length(epochIn)+1));
                segLength=get(epochOut,'segLength');
                set(epochOut,'segLength',cat(2,segLength,length(seg)));
            end

        end





        function S=nl2seg(N, varName)
            % Convert an nldat object to a segdat using nans as segment separators;
            % S is empty if there is no valid data in N

            [nSamp,nChan,nReal]=size(N);
            if nChan>1,
                d=N.dataSet;
                iNAN=find(isnan(d(:,1)));
                for iChan=2:nChan
                    jNAN=find(isnan(d(:,iChan)));
                    if iNAN~=jNAN
                        warning('Channels have different segmentation ');
                    end
                end
            end
            %% Generate data set for segment analysis
            dN=domain(N);
            dStart=min(dN);
            dEnd=max(dN);
            dIncr=N.domainIncr;
            nSamp=1+round((dEnd-dStart)/N.domainIncr);
            dataSet=nan(nSamp,nChan);
            idx=idx4domain(dStart,N.domainIncr, dN);
            dataSet(idx,:)=N.dataSet;
            nDomain=[dStart:dIncr:dEnd];
            % Find start and end points of data
            segCnt=0;
            segStart=1;
            nLen=length(dataSet);
            % Find start and end of segments
            c=categorical;
            c(1:nLen)='good';
            iNan=find(isnan(dataSet(:,1)));
            c(iNan)='nan';
            e=eseq(c);
            ne=length(e);
            seqCnt=0;
            segIno={};
            onsetPointer=0;
            newDataSet=[];
            %nDomain=domain(N);
            newComment={};
            for ie=1:ne,
                if e(ie).type=='good'
                    segCnt=segCnt+1;
                    domainStart(segCnt)=nDomain(e(ie).startIdx);
                    onsetPointer(segCnt)=length(newDataSet)+1;
                    segLength(segCnt)=e(ie).nSamp;
                    segInfo{segCnt}=[ varName num2str(segCnt)];
                    chanNames=N.chanNames;
                    t=chanNames{1};
                    newDataSet=cat(1,newDataSet, dataSet(e(ie).startIdx:e(ie).endIdx,:));
                end
            end

            S=segdat;


            if length(newDataSet)>0
                set(S,'chanNames',N.chanNames, 'chanUnits',N.chanUnits, 'domainIncr',N.domainIncr, ...
                    'domainStart',domainStart,'domainValues',nan, 'dataSet', newDataSet, ...
                    'dataSize', size(dataSet),'comment',N.comment, ...
                    'onsetPointer', onsetPointer,'segLength',segLength, 'segInfo',segInfo);
            else
                set(S,'dataSet',[]);
            end
        end

        function S = randSeg(Z, minSegLen, maxSegLen, maxSegSpace)
            % randSeg - generate a random segdat from a nldat object
            % S = randSeg(Z, minSegLen, maxSegLen, maxSegSpace)
            % Z - nldat object
            % minSegLen - minimum segment length
            % maxSegLen - maximum sewgment length
            % maxSegSpace - maximum spacing between segments
            % S - segdat object sampled randomly from Z
            nMax=length(Z);
            iEnd=0;
            continueFlag=true;
            S=segdat;
            i=1;
            while continueFlag
                segLen=randi([minSegLen maxSegLen]);
                segOffSet=randi([1 maxSegSpace]);
                iStart=iEnd + segOffSet;
                iEnd=iStart+segLen-1;
                if iEnd>nMax,
                    iEnd=nMax;
                    continueFlag=false;
                end
                sTemp=segdat(Z(iStart:iEnd,:));
                if i==1
                    S=sTemp;
                else
                    S=cat(1,S,sTemp);
                end
                i=i+1;

            end
        end



    end
end



