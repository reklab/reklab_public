classdef ssm<nltop
    %State-space model for NLID toolbox
    properties
        parameterSet
    end
    methods
        function S = ssm (a,varargin)
            %SSM specific parameters
            S.parameterSet=param('paramName','A','paramDefault',[], ...
                'paramHelp','State-Space matrix A', ...
                'paramType','number');
            S.parameterSet(2)=param('paramName','B','paramDefault',[], ...
                'paramHelp','State-Space matrix A', ...
                'paramType','number');
            S.parameterSet(3)=param('paramName','C','paramDefault',[], ...
                'paramHelp','State-Space matrix A', ...
                'paramType','number');
            S.parameterSet(4)=param('paramName','D','paramDefault',[], ...
                'paramHelp','State-Space matrix A', ...
                'paramType','number');
            S.parameterSet(5)=param('paramName','idMethod','paramDefault','PI', ...
                'paramHelp','Subspace id Method (PI/PO) ', ...
                'paramType','select','paramLimits', {'PI' 'PO' } );
            S.parameterSet(6)=param('paramName','order','paramDefault',2, ...
                'paramType','number',  ...
                'paramHelp','System orde');

            S.parameterSet(7)=param('paramName','orderSelect','paramDefault','manual', ...
                'paramType','select', 'paramLimits',{'manual' 'largest-gap' 'preset'}, ...
                'paramHelp','Methods for order seletion');
            S.parameterSet(8)=param('paramName','hankleSize','paramDefault',20, ...
                'paramHelp','Size of hankle matrix (> system order)', ...
                'paramType','number');
            S.parameterSet(9)=param('paramName','domainIncr','paramDefault',0.001, ...
                'paramHelp','Sampling Time', ...
                'paramType','number');
            S.parameterSet(10)=param('paramName','nDelayInput','paramDefault',0, ...
                'paramHelp','Input delay in samples', ...
                'paramType','number');
            S.parameterSet(11)=param('paramName','displayFlag','paramDefault',0,...
                'paramHelp','display','paramLimits', {0,1});
            set(S,'comment','State-Space Model');
            if nargin==0;
                return
            elseif nargin==1,
                S=nlmkobj(S,a);
            elseif isa(a,'ssm')
                S = nlmkobj(a,varargin{:});
            else
                S = nlmkobj(S,a,varargin{:});
            end
            
        end
        function I = irf(S)
            % Convert an SSM to IRF
            % Note that in this version, ssm does not represent non-causal system
            % So, conversion to double sided IRF is not possible
            I = irf;
            matrix_a = get(S,'A');
            matrix_b = get(S,'B');
            matrix_c = get(S,'C');
            matrix_d = get(S,'D');
            delay = get(S,'nDelayInput');
            ts = get(S,'domainIncr');
            sys_ss = ss(matrix_a,matrix_b,matrix_c,matrix_d,ts,'InputDelay',delay);
            [~,lag_max] = impulse(sys_ss);
            lag_max = lag_max(end);
            time = 0 : ts : lag_max;
            i = impulse(sys_ss,time);
            set(I,'dataSet',i,'domainIncr',ts);
        end
        function F = fresp(S)
            % Convert an SSM to fresp
            F = fresp;
            nFFT = 500;
            matrix_a = get(S,'A');
            matrix_b = get(S,'B');
            matrix_c = get(S,'C');
            matrix_d = get(S,'D');
            delay = get(S,'nDelayInput');
            ts = get(S,'domainIncr');
            sys_ss = ss(matrix_a,matrix_b,matrix_c,matrix_d,ts,'InputDelay',delay);
            frequency_axis = 2*pi*linspace(0,1/2/ts,nFFT);
            freq_incr = (frequency_axis(2)-frequency_axis(1))/2/pi;
            [mag,phase] = bode(sys_ss,frequency_axis);
            mag = shiftdim(mag,2);
            phase = shiftdim(phase,2);
            z = mag.*exp(phase*pi/180*sqrt(-1));
            set(F,'dataSet',z,'domainStart',0, 'domainIncr',freq_incr);
        end
        function plot(S)
            I = irf;
            matrix_a = get(S,'A');
            matrix_b = get(S,'B');
            matrix_c = get(S,'C');
            matrix_d = get(S,'D');
            delay = get(S,'nDelayInput');
            ts = get(S,'domainIncr');
            sys_ss = ss(matrix_a,matrix_b,matrix_c,matrix_d,ts,'InputDelay',delay);
            [~,lag_max] = impulse(sys_ss);
            lag_max = lag_max(end);
            time = 0 : ts : lag_max;
            i = impulse(sys_ss,time);
            set(I,'dataSet',i,'domainIncr',ts,'comment',['IRF Representation of ',get(S,'comment')]);
            plot(I);
            
        end
        function out = nlsim(S,z)
            [~,nchan,~]=size(z);
            ztype = class(z);
            switch ztype
                case 'segdat'
                    disp('ssm.nlsim requires input and output for accurate simulation of segdat data.');
                    
                    if nchan == 2
                        out = nlsim_short_segment (S.parameterSet,z)
                    end
                        case 'nldat'
                            matrix_a = get(S,'A');
                            matrix_b = get(S,'B');
                            matrix_c = get(S,'C');
                            matrix_d = get(S,'D');
                            ts = get(S,'domainIncr');
                            delay = get(S,'nDelayInput');
                            input = get(z,'dataSet');
                            inputd = del(input,delay);
                            out = dlsim(matrix_a,matrix_b,matrix_c,matrix_d,inputd);
                            out = nldat(out,'domainIncr',ts,'domainStart',z.domainStart);
                            otherwise
                                error('The input type not supported.')
                    end
            end
            function S = nlident (S, z,  varargin)
                %Identify a state-space object from input-output data
                %
                if nargin < 2,
                    disp('NLIDtakes three inputs for SSM objects: irf, fresp, nldat' );
                elseif nargin > 2,
                    set(S,varargin{:});
                end
                assign(S.parameterSet)
                set(S,'domainIncr',get(z,'domainIncr'));
                ztype = class(z);
                switch ztype
                    case 'fresp'
                        warning('This option is not implemented!');
                        S = ssm;
                    case 'irf'
                        warning('This option is not implemented!');
                        S = ssm;
                    case 'segdat'
                        S = short_segment(z,S.parameterSet);
                    case 'nldat'
                        in = z(:,1);
                        out = z(:,2);
                        in = delay(in,nDelayInput);
                        x = get(in,'dataSet');
                        y = get(out,'dataSet');
                        switch upper(idMethod)
                            case 'PI'
                                % extract the column space using the i past input moesp algorithm
                                %which can be used for the output error identification problem. 
                                [Sn,R] = dordpi(x,y,hankleSize);
                                if ~strcmp(orderSelect,'preset'),
                                    order = orderselect(Sn,orderSelect);
                                end
                                [A,C]  =  destac(R,order);
                                [B,D] =  destbd(x,y,A,C);
                                if any(isnan(B)) | any(isnan(D))
                                set(S,'comment','ID failed', 'order',-1) % Return an empty value if ID failed. 
                                    return;
                                end
                                set(S,'A',A,'B',B,'C',C,'D',D,'order',order);
                                if get(S,'displayFlag') == 1
                                    predicted_output = nlsim(S,in);
                                    figure;
                                    plot(cat(2,out,predicted_output),'plotmode','super');
                                end
                            case 'PO'
                                % extract the column space usibg the past output moesp algorithm  which can be 
                                % used for the innovations model identification  problem. 
                               [Sn,R] = dordpo(x,y,hankleSize);
                                if ~strcmp(orderSelect,'preset')
                                  order = orderselect(Sn,orderSelect);
                                end
                                [A,C]  =  destac(R,order);
                                [B,D] =  destbd(x,y,A,C);
                                if isempty (B)
                                    B=0;
                                end
                                if isempty (D)
                                    D=0;
                                end
                                set(S,'A',A,'B',B,'C',C,'D',D,'order',order);
                                if get(S,'displayFlag') == 1
                                    predicted_output = nlsim(S,in);
                                    figure;
                                    plot(cat(2,out,predicted_output),'plotmode','super');
                                end
                            otherwise
                                error('Identification method must be PI or PO');
                        end
                    otherwise
                        error (['ssm nlident - does not support data type: ' ztype]);
                end
            end
        end
    end
    
