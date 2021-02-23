classdef eseq
    % UESEQ - event sequence oject
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.
    
    % Public, tunable properties
    properties
        domainStart=0;
        domainIncr=1;
        startIdx=[];
        endIdx=[];
        nSamp=[];
        type = categorical
        
    end
    
    
    
    
    methods
        function e = eseq(a)
            if nargin==0
                return
            elseif isa(a,'categorical')
                e=eseq.cseq2eseq(a);
            else
                error ('Invlaid input type');
            end
            
            
            
        end
        
        function [cseq, domainVal] = cseq (eseq)
            % eseq2cseq - converts an event sequence to to a categorical sequence
            
            cseq=categorical;
            idxMin=min([eseq.startIdx]);
            idxMax=max([eseq.endIdx]);
            cseq=categorical;
            
            d=[];
            for i=1:length(eseq)
                nCur=eseq(i).endIdx-eseq(i).startIdx+1;
                curC=categorical;
                curC(1:nCur,1)=eseq(i).type;
                cseq=cat(1,cseq,curC);
                curD=domain(eseq(i));
                d=cat(1,d,[curD{:}]');
            end
            if nargout==2,
                domainVal=d;
            end
        end
        
        function d=domain(eseq)
            % return a cell array of domain values for an event sequence
            d={};
            nEvent=length(eseq);
            for iEvent=1:nEvent,
                eCur=eseq(iEvent);
                x=[eCur.startIdx:eCur.endIdx]-1;
                curDomain=eCur.domainStart+x*eCur.domainIncr;
                d{iEvent}=curDomain;
            end
        end
        
        function eInter = intersect (e1,e2);
            % return events where e1 and e2  are of the same type and intersect
            % assumes that e1 and e2 are in increaeing time and msut have
            % same domainStart and domainIcrc=
            eInter=eseq;
            n1=length(e1);
            n2=length(e2);
            d1=domain(e1);
            d2=domain(e2);
            iInter=0;
            for i1=1:n1,
                e1Cur=e1(i1);
                e1CurStart=e1Cur.startIdx
                e1CurEnd=e1Cur.endIdx;
                for i2=1:n2,
                    e2Cur=e2(i2);
                    e2CurStart=e2Cur.startIdx;
                    e2CurEnd=e2Cur.endIdx;
                    if e1CurEnd<e2CurStart | e2CurEnd<e1CurStart
                        continue
                    else
                        if e1Cur.type==e2Cur.type
                        iInter=iInter+1;
                        eInter(iInter,1).startIdx=max(e1CurStart, e2CurStart);
                        eInter(iInter,1).endIdx=min(e1CurEnd,e2CurEnd);
                        eInter(iInter,1).type=e1Cur.type;
                        end
                    end
                end
            end
        end
             function h=line(e)
                [c,d]=cseq(e);
                line(d,c);
                
             end
            
            
            function plot(e)
                [c,d]=cseq(e);
                plot (d,c,'o');
                
            end
        end

        
        methods (Static)
            function event =cseq2eseq(cSeq)
                %catSignal2Events - convert a categorical sequence  to event sequence
                % Input: cSig - categorical signal
                % output
                % event - structure with fields
                % event.type
                % event.startIdx
                % event.endIdx,
                % event.length - event length in samples
                if iscell(cSeq),
                    n=length(cSeq);
                    for i=1:n,
                        e=cseq2eseq(cSeq{i}, i);
                        event{i,1}=e;
                    end
                else
                    
                    
                    if nargin==1,
                        caseNum=-1;
                    end
                    
                    nSamp=length(cSeq);
                    x=double(cSeq);
                    diffX=diff(x);
                    iChange=find(diffX~=0);
                    nEvent=length(iChange);
                    event=eseq;
                    
                    if nEvent==0,
                        event(1).startIdx=1;
                        event(1).endIdx=length(cSeq);
                        event(1).type=cSeq(1) ;
                        event(1).nSamp=event(1).endIdx-event(1).startIdx+1;
                        return
                    end
                    
                    
                    event(1).startIdx=1;
                    event(1).endIdx=iChange(1);
                    event(1).type=cSeq(1) ;
                    event(1).nSamp=iChange(1);
                    
                    
                    for iEvent=2:nEvent
                        event(iEvent).startIdx=event(iEvent-1).endIdx+1;
                        event(iEvent).endIdx=iChange(iEvent);
                        event(iEvent).type=cSeq(iChange(iEvent));
                        event(iEvent).nSamp=event(iEvent).endIdx-event(iEvent).startIdx+1;
                    end
                    
                    event(nEvent+1).startIdx=iChange(nEvent)+1;
                    event(nEvent+1).endIdx=length(cSeq);
                    event(nEvent+1).type=cSeq(end);
                    event(nEvent+1).nSamp=event(nEvent+1).endIdx-event(nEvent+1).startIdx+1;
                    event=event';
                    
                    
                end
            end
        end
    end
    
