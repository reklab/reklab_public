classdef segdat<nldat
    %Segmented nldat object
    properties
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
                S=nl2seg(a);
            else
                set (S,{a varargin{:}});
            end
            
        end
        
        
        
        function Z = cat (DIM, varargin)
            if nargin <3.
                error('cat must have at least three input arguments');
            end
            Z = varargin{1};
            [nSamp,nChan,nReal]=size(Z);
            z = double(varargin{1});
            onsetpointer = double(get(Z,'onsetPointer'));
            seglength = double(get(Z,'segLength'));
            domainStart=get(Z,'domainStart');
            Z.comment = [ 'cat' int2str(DIM) ' ' inputname(2)  ];
            
            for i = 2 : nargin - 1
                curZ=varargin{i};
                [nSampCur, nChanCur, nRealCur]=size(curZ);
                if DIM==1,
                    Z=segCat(Z,varargin{i});
                elseif DIM==2
                    if nSamp ~= nSampCur
                        error ('Number of samples must be equal if adding channels');
                    end
                    Z.dataSet=cat(2,Z.dataSet,curZ.dataSet);
                    Z.chanNames=cat (1,Z.chanNames,curZ.chanNames);
                else
                    error ([' DIM = ' numstr(DIM) ' not yet supported']);
                end
                Z.comment = [Z.comment ',' inputname(i+1)];
                
            end
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
                    plot(sSeg,'lineColor',colors(i,:))
                end
             
            else
                plot(nldat(S))
            end
        end
        
        function Z = segCat (Z1, Z2)
            % Concatonate segdat objects - only for dim 1.
            % In case of overlap the output is given by Z2
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
            zd(ptr1,:)=Z1.dataSet(:,:);
            ptr2=idx4domain(dMin,domainIncr,d2);
            iOverlap=intersect(ptr1,ptr2);
            comment=['segcat(' name1 ',' name2 ')'];
            if ~isempty (iOverlap),
                overlapStart=min(iOverlap);
                overlapEnd=max(iOverlap);
                comment=['segcat:' name2 ' overlaps ' name1 ' from ' num2str(overlapStart) ...
                    ':' num2str(overlapEnd)];
            end
            zd(ptr2,:)=Z2.dataSet(:,:);
            set (n1,'domainStart',dMin,'dataSet',zd, 'comment' ,comment);
            Z=segdat(n1);        
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
                subplot (nChan,1,iChan);
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
                'domainValues',nan,'comment','nldat(segdat)');
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
            for k = 1: nchan
                for i = 1 : length(onsetpointer)
                    for j = 1 : length(onsetpointer)
                        if i == j
                            break;
                        elseif (onsetpointer(j)<=onsetpointer(i))&&(onsetpointer(i)<=endpointer(j))
                            warning(['In channel ',num2str(k),', segment ',num2str(i),' & segment',num2str(j),' overlap.'])
                        end
                    end
                end
                
            end
            
        end
        
    end
    
end

function S=nl2seg(N)
% Convert an nldat object to a segdat using nans as segment separators;
[nSamp,nChan,nReal]=size(N);
if nChan>1,
    error('segdat does not yet support mulitple channels');
end
dataSet=N.dataSet;
% Find start and end points of data
segCnt=0;
segStart=1;
nLen=length(dataSet);
% Find start and end of segments
for iCur=1:nLen-1,
    curVal=dataSet(iCur,1);
    nextVal=dataSet(iCur+1,1);
    if iCur==1 & ~isnan(curVal) % Start of first section
        segCnt=1;
        segStart(segCnt)=iCur;
    elseif isnan(curVal) & ~isnan(nextVal) % start of a section
        segCnt=segCnt+1;
        segStart(segCnt)=iCur+1;
    elseif ~isnan(curVal) & isnan(nextVal) % End of a section
        segEnd(segCnt)=iCur;
    elseif iCur==nLen-1;
        segEnd(segCnt)=iCur+1;
    end
end
nDomain=domain(N);
domainStart=nDomain(segStart);
segLength=segEnd-segStart+1;
onsetPointer(1)=1;
for i=2:segCnt,
    onsetPointer(i)=onsetPointer(i-1)+segLength(i-1);
end
S=segdat;
i=find(isnan(dataSet));
dataSet(i,:)=[];
set(S,'chanNames',N.chanNames, 'chanUnits',N.chanUnits, 'domainIncr',N.domainIncr, ...
    'domainStart',domainStart,'domainValues',nan, 'dataSet', dataSet, ...
    'dataSize', size(dataSet),'comment',N.comment, ...
    'onsetPointer', onsetPointer,'segLength',segLength);
end

