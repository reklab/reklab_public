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
                    S.dataSet = a(:);
                    if nargin > 1,
                        set (S,varargin{:});
                    end
            elseif isa(a,'segdat')
                S = nlmkobj(a,varargin{:});
            else
                S = nlmkobj(S,a,varargin{:});
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
            Z.comment = [ 'cat' int2str(DIM) ' ' inputname(2)  ];
            
            for i = 2 : nargin - 1
                y = varargin{i};
                z = cat(DIM, z , double(y));
                onsetpointer = cat(DIM,onsetpointer,get(y,'onsetPointer'));
                seglength = cat(DIM,seglength,get(y,'segLength'));
                Z.comment = [Z.comment ',' inputname(i+1)];
                switch DIM
                    case 2
                        Z.chanNames=cat(2,Z.chanNames,y.chanNames);
                end
            end
            set(Z,'dataSet',z,'onsetPointer',onsetpointer,'segLength',seglength);
        end
        function plot(S)
            if (size(S,2)==1) && ((size(S,3)==1))
                S_nldat = nldat(get(S,'dataSet'),'domainIncr',get(S,'domainIncr'),'chanNames'...
                ,get(S,'chanNames'),'chanUnits',get(S,'chanUnits'),'domainName',get(S,'domainName')...
                ,'comment',get(S,'comment'));
                plot(S_nldat)
                hold on
                numSegment = length(get(S,'segLength'));
                onsetPointer = get(S,'onsetPointer');
                segLength = get(S,'segLength');
                for i = 1 : numSegment
                    S_nldat_Segment = S_nldat(onsetPointer(i):onsetPointer(i)+segLength(i)-1,:);
                    plot(S_nldat_Segment,'line_color','r')
                end
                legend('Original Record','Segmented Data')
            else
                plot(nldat(S))
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
            dataset = get(data,'dataSet');
            d = zeros(sum(seglength(:,1)),nchan);
            for i = 1: nchan
                pointer = 1;
                for j = 1 : size(onsetpointer,1)
                     onsetpointer_new(j,i) = pointer;
                     seglength_new(j,i) = length(dataset(onsetpointer(j,i):endpointer(j,i),i));
                     d(onsetpointer_new(j,i):onsetpointer_new(j,i)+seglength_new(j,i)-1,i) = dataset(onsetpointer(j,i):endpointer(j,i),i);
                     pointer = pointer + seglength_new(j,i);
                end
            end
            out = nldat(d,'chanNames',data.chanNames,'chanUnits',data.chanUnits,'domainIncr',data.domainIncr,'domainName',data.domainName,'domainValues',data.domainValues,'dataSize',data.dataSize,'comment',data.comment);
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
                error('At least one semgnet does not have a correct length')
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

