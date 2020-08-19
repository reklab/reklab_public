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
                set(S,'onsetPointer',ones(1,nChan));
                set(S,'segLength',nSamp*ones(1,nChan));
                
                if nargin > 1,
                    set (S,varargin{:});
                end
            elseif isa(a,'segdat')
                S = nlmkobj(a,varargin{:});
            elseif isa(a,'nldat');
                fields=fieldnames(a);
                nField=length(fields);
                for i=1:nField,
                    S.(fields{i})=a.(fields{i});
                end
                [nSamp,nChan]=size(a);
                set(S,'onsetPointer',1);
                set (S,'segLength', nSamp);
                if nargin > 1,
                    set (S,varargin{:});
                end
            end
            
        end
        function Z = cat (DIM, varargin)
            if nargin <3.
                error('cat must have at least three input arguments');
            end
            Z = varargin{1};
            z = double(varargin{1});
            onsetpointer = double(get(Z,'onsetPointer'));
            seglength = double(get(Z,'segLength'));
            domainStart=get(Z,'domainStart'); 
            Z.comment = [ 'cat' int2str(DIM) ' ' inputname(2)  ];
            
            for i = 2 : nargin - 1
                y = varargin{i};
                z = cat(DIM, z , double(y));
                onsetpointer = cat(DIM,onsetpointer,newOnset);
                seglength = cat(DIM,seglength,get(y,'segLength'));
                Z.comment = [Z.comment ',' inputname(i+1)];
                switch DIM
                    case 1
                        domainStart=cat(1,domainStart,y.domainStart);
                    case 2
                        Z.chanNames=cat(2,Z.chanNames,y.chanNames);
                end
            end
            set(Z,'dataSet',z,'onsetPointer',onsetpointer,'segLength',seglength, ...
                'domainStart',domainStart);
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
                legend('Original Record','Segmented Data')
            else
                plot(nldat(S))
            end
        end
        
        function Z = segCat (Zin, varargin)
            % Concatonate segdat objects - only for dim 1.
            % only works for nonoverlapping segments
            % error checking very weak
            if isa(Zin,'nldat');
                Z=segdat(Zin);
            else
                Z=Zin;
            end
                
            nSeg=nargin-1;
            for i=1:nSeg
                V=varargin{i};
                if isa(V,'nldat'),
                   V=segdat(V);
                end
                Z.dataSet=cat(1,Z.dataSet, V.dataSet);
                Z.domainStart=cat(1,Z.domainStart,V.domainStart);
                curSegLen=get(Z,'segLength');
                newSegLen=get(V,'segLength');
                updateSegLen=cat(1,curSegLen, newSegLen);
                curOnset=get(Z,'onsetPointer');
                newOnset=sum(curSegLen,1)+curOnset;
                updateOnset=cat(1,curOnset,newOnset); 
                set(Z,'onsetPointer',updateOnset,'segLength',updateSegLen);
            end
        end
                
        
        function n = segCount(S)
            % reutrns number of segments in a segdat object
            onsetPointer=get(S,'onsetPointer');
            n=length(onsetPointer);
        end
        
        function N = segGet(S,segNum)
            % Return segment of segdat object as a nldat object
            nSeg=segCount(S);
            if (segNum>nSeg)
                error('too many segments')
            end
            nTmp=nldat(S);
            segLen=get(S,'segLength');
            
            segStartPointer=[ 0 ; cumsum(segLen(:,1))]+1;
            sStart=segStartPointer(segNum);
            sEnd=segStartPointer(segNum+1)-1;            
            N=nTmp(sStart:sEnd,:);
            N.comment=(['Segment ' num2str(segNum) '  of segdat object']); 
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
                    d_temp = decimate(dataset(onsetpointer(j,i):endpointer(j,i),i),decimation_ratio);
                    onsetpointer_new(j,i) = pointer;
                    seglength_new(j,i) = length(d_temp);
                    d(onsetpointer_new(j,i):onsetpointer_new(j,i)+seglength_new(j,i)-1,i) = d_temp;
                    pointer = pointer + seglength_new(j,i);
                end
            end
            set(out,'dataSet',d,'onsetPointer',onsetpointer_new,'segLength',seglength_new);
        end
        function out = ddt(data)
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
            dataset = get(data,'dataSet');
            d = zeros(sum(seglength(:,1)),nchan);
            for i = 1: nchan
                pointer = 1;
                for j = 1 : size(onsetpointer,1)
                    d_temp = dataset(onsetpointer(j,i):endpointer(j,i),i);
                    d_temp = nldat(d_temp,'domainIncr',ts);
                    d_temp = ddt(d_temp);
                    d_temp = d_temp.dataSet;
                    onsetpointer_new(j,i) = pointer;
                    seglength_new(j,i) = length(d_temp);
                    d(onsetpointer_new(j,i):onsetpointer_new(j,i)+seglength_new(j,i)-1,i) = d_temp;
                    pointer = pointer + seglength_new(j,i);
                end
            end
            set(out,'dataSet',d,'onsetPointer',onsetpointer_new,'segLength',seglength_new);
        end
        function out = nldat(data)
            % returns concatonated segments as a nldat object
            
            errorcheck(data);
            [~,nchan,~] =size(data);
            onsetpointer = get(data,'onsetPointer');
            seglength = get(data,'segLength');
            if size(seglength,2)==1
            elseif  any(~(mean(seglength')==seglength(:,1)'))
                error('All channels in each  data segment must have equal segment lengths')
            end
            nSeg=segCount(data);
            [nSamp,nChan,nReal]=size(data);
            domainIncr=data.domainIncr; 
            newDomain=[];; 
            j=0;
            if nSeg==1,
                newDomain=nan;
            else
            for iSeg=1:nSeg,
                iStart=onsetpointer(iSeg);
                curSegLen=seglength(iSeg);
                domainStart=data.domainStart(iSeg);
                xTemp=domainStart+([1:curSegLen]-1)*domainIncr;
                newDomain=cat(1,newDomain,xTemp(:));
            end
            end
            dataset = get(data,'dataSet');
           out = nldat(dataset,'chanNames',data.chanNames,'chanUnits',data.chanUnits, ...
                'domainIncr',data.domainIncr,'domainName',data.domainName,'domainStart',data.domainStart, ...
                'domainValues',newDomain, ...
                'dataSize',data.dataSize,'comment',data.comment);
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
            if ~(size(onsetpointer,1) == size(seglength,1))
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
                for i = 1 : length(onsetpointer(:,k))
                    for j = 1 : length(onsetpointer(:,k))
                        if i == j
                            break;
                        elseif (onsetpointer(j,k)<=onsetpointer(i,k))&&(onsetpointer(i,k)<=endpointer(j,k))
                            warning(['In channel ',num2str(k),', segment ',num2str(i),' & segment',num2str(j),' overlap.'])
                        end
                    end
                end
                
            end
            
        end
    end
end

