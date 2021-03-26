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
                S.dataSet = a;
                [nSamp,nChan]=size(a);
                if nChan>1,
                    error('segdat: Multiple channels not yet supported');
                end
                set(S,'onsetPointer',1);
                set(S,'segLength',nSamp);
                if nargin > 1,
                    set (S,varargin{:});
                end
            elseif isa(a,'segdat')
                S = nlmkobj(a,varargin{:});
            elseif isa(a,'nldat')
                S=nl2seg(a,inputname(1));
            else
                set (S,{a varargin{:}});
            end
            
        end
        
        function N=seg4domain (S, domainVal)
            % Return segment mumber associated with a domain Value
            dStart=S.domainStart;
            dEnd=domainEnd(S);
            N=find(domainVal>=dStart & domainVal<=dEnd);
        end
        
        
        function sCat = cat (S1,S2)
            % Concatonate two segdat objects
            domainIncr=S1.domainIncr;
            d1=domain(S1);
            s1Data=double(S1);
            d2=domain(S2);
            s2Data=double(S2);
            catVector =categorical;
            d2=domain(S2);
            domainStart=min(min(d1),min(d2));
            domainEnd=max(max(d1),max(d2));
            
            idxMax=idx4domain(domainStart, domainIncr, domainEnd);
            domainVector=domainStart:domainIncr:domainEnd;
            valueVector(:,1)=domainVector*0;
            catVector(1:idxMax)='0';
            idx1=idx4domain(domainStart, domainIncr, d1);
            valueVector(idx1,1)=s1Data;
            catVector(idx1)='1';
            idx2=idx4domain(domainStart, domainIncr, d2);
            valueVector(idx2,1)=s2Data;
            catVector(idx2)='2';
            idxIntersect=intersect(idx1,idx2);
            if ~isempty(idxIntersect)
                catVector(idxIntersect)='3';
            end
            e=eseq(catVector);
            
            eDomain=domain(e,domainStart,domainIncr);
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
                    dataSet=cat(1, dataSet, valueVector(iStart:iEnd,1));
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
        end
        
        
        function plot(S)
            colors=colororder;
            if (size(S,2)==1) && ((size(S,3)==1))
                S_nldat=nldat(S);
                % plot(S_nldat)
                hold on
                numSegment = segCount(S);
                onsetPointer = get(S,'onsetPointer');
                segLength = get(S,'segLength');
                for i = 1 : numSegment
                    sSeg=segGet(S,i);
                    plot(sSeg)
                end
                
            else
                plot(nldat(S))
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
        
        
        
        function Z = segIntersect (Z1, Z2)
            % Determine intersection between two segdat objects
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
            zd(ptr1,:)=1;
            ptr2=idx4domain(dMin,domainIncr,d2);
            zd(ptr2,:)=2;
            iOverlap=intersect(ptr1,ptr2);
            comment=['segcat(' name1 ',' name2 ')'];
            if ~isempty (iOverlap),
                zd(iOverlap,:)=3;
            end
            Z=nldat (n1,'domainStart',dMin,'dataSet',zd );
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
            dataSet=S.dataSet;
            segLen=get(S,'segLength');
            onsetPointer=get(S,'onsetPointer');
            sStart=onsetPointer(segNum);
            sEnd=onsetPointer(segNum)+segLen(segNum)-1;
            N.dataSet=dataSet(sStart:sEnd,:);
            set(N, 'chanNames', S.chanNames,'chanUnits',S.chanNames, 'domainIncr',S.domainIncr, ...
                'domainName',S.domainName);
            N.comment=(['sSegment ' num2str(segNum) '  of segdat object']);
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
            if size(seglength,2)==1
            elseif  any(~(mean(seglength')==seglength(:,1)'))
                error('All data channels must have equal segment lengths')
            end
            onsetpointer_new = zeros(size(onsetpointer));
            seglength_new = zeros(size(onsetpointer));
            endpointer = onsetpointer + seglength - 1;
            ts = get(data,'domainIncr');
            out = data;
            set(out,'domainIncr',ts*decimation_ratio);
            dataset = get(data,'dataSet');
            d = zeros(ceil(sum(seglength(:,1))/decimation_ratio),nchan);
            for i = 1: nchan
                pointer = 1;
                for j = 1 : size(onsetpointer,1)
                    d_temp = decimate(dataset(onsetpointer(j):endpointer(j),i),decimation_ratio);
                    onsetpointer_new(j) = pointer;
                    seglength_new(j) = length(d_temp);
                    d(onsetpointer_new(j):onsetpointer_new(j)+seglength_new(j)-1) = d_temp;
                    pointer = pointer + seglength_new(j);
                end
            end
            set(out,'dataSet',d,'onsetPointer',onsetpointer_new,'segLength',seglength_new);
        end
        
        function out = ddt(data)
            errorcheck(data);
            nSeg=segCount(data);
            for iSeg=1:nSeg
                xTmp=segGet(data,iSeg);
                xdt=ddt(xTmp);
                if iSeg==1,
                    out=segdat(xdt);
                else
                    out=segCat(out,xdt);
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
            data1 = nldat(data1);
            data2 = nldat(data2);
            v = vaf(data1,data2);
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
        
    end
    
end

function S=nl2seg(N, varName)
% Convert an nldat object to a segdat using nans as segment separators;

[nSamp,nChan,nReal]=size(N);
if nChan>1,
    error('segdat does not yet support mulitple channels');
end
dataSet=N.dataSet;
% Find start and end points of datas
segCnt=0;
segStart=1;
nLen=length(dataSet);
% Find start and end of segments
c=categorical;
c(1:nLen)='good';
iNan=find(isnan(dataSet));
c(iNan)='nan';
e=eseq(c);
ne=length(e);
seqCnt=0;
segIno={};
onsetPointer=0;
newDataSet=[];
nDomain=domain(N);
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
        
        newDataSet=cat(1,newDataSet, dataSet(e(ie).startIdx:e(ie).endIdx));
    end
end

S=segdat;
i=find(isnan(dataSet));
dataSet(i,:)=[];
set(S,'chanNames',N.chanNames, 'chanUnits',N.chanUnits, 'domainIncr',N.domainIncr, ...
    'domainStart',domainStart,'domainValues',nan, 'dataSet', newDataSet, ...
    'dataSize', size(dataSet),'comment',N.comment, ...
    'onsetPointer', onsetPointer,'segLength',segLength, 'segInfo',segInfo);
end

