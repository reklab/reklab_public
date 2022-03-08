classdef eseq
    % eseq - class for manipulation of event sequences.
    %
    %
    
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
        function e = eseq(a, domainStart,domainIncr, eventCat)
            if nargin==0
                return
            end
            if nargin<3
                domainIncr=1;
            end
            if nargin<2
                domainStart=0;
            end
            
            if isa(a,'categorical')
                e=eseq.cseq2eseq(a, domainStart, domainIncr);
            elseif isa(a,'double')
                e=eseq.idx2eseq(a, domainStart,domainIncr, eventCat);
                
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
        
        
        function d=domain(e)
            %  d=domain(eseq) return a cell array of domain values for an event sequence
            d={};
            nEvent=length(e);
            for iEvent=1:nEvent,
                eCur=e(iEvent);
                x=[eCur.startIdx:eCur.endIdx]-1;
                curDomain=eCur.domainStart+x*eCur.domainIncr;
                d{iEvent}=curDomain;
            end
        end
        
        function  eOut=events4epoch ( e, epochStart, epochEnd, overlapFlag, minOverlapLen)
            startIdx=[e.startIdx];
            endIdx=[e.endIdx];
            domainStart=[e.domainStart];
            domainIncr=[e.domainIncr];
            startTime=domainStart +(startIdx-1).*domainIncr;
            endTime=domainStart +(endIdx-1).*domainIncr;
            if ~overlapFlag
            i= find(startTime>= epochStart & endTime <= epochEnd);
            else
                i= (startTime>= epochStart & endTime <= epochEnd) | ...
                    (startTime<epochStart & (endTime-epochStart)>=minOverlapLen) | ...
                    (endTime>epochEnd & (epochEnd-startTime)>=minOverlapLen) | ...
                    (startTime < epochStart & endTime > epochEnd);                
            end
            eOut=e(i);
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
                e1CurStart=e1Cur.startIdx;
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
                            eInter(iInter,1).domainStart=e1.domainStart;
                            eInter(iInter,1).domainIncr=e1.domainIncr;
                            eInter(iInter,1).startIdx=max(e1CurStart, e2CurStart);
                            eInter(iInter,1).endIdx=min(e1CurEnd,e2CurEnd);
                            eInter(iInter,1).nSamp=  eInter(iInter,1).endIdx - eInter(iInter,1).startIdx +1
                            eInter(iInter,1).type=e1Cur.type;
                        end
                    end
                end
            end
        end
        
        function nOut= nldat4eseq ( E, N)
           %  N = nldat4eseq ( E, N)
           % E - event sequence arrau
           % nldat or segdat object
           % Returns a cell array of nldat objects one for each event in E.
           % Each nldat object containts the data from N correspond to the
           % event.
           if isa(N,'segdat'),
               N=nldat(N);
           end
           nOut={};
           eDomain=domain(E); 
           nEvent=length(E);
           for iEvent=1:nEvent,
               curDomain=eDomain{iEvent};
               idx=idx4domain(N.domainStart, N.domainIncr, curDomain);
               curOut= N(idx,:);
               nOut{iEvent}=curOut;              
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
        function event =cseq2eseq(cSeq, domainStart, domainIncr)
            %catSignal2Events - convert a categorical sequence  to event sequence
            % Input: cSig - categorical signal
            % output
            % event - structure with fields
            % event.type
            % event.startIdx
            % event.endIdx,
            % event.length - event length in samples
            if nargin<3
                domainIncr=1;
            end
            if nargin<2
                domainStart=0;
            end
            
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
                    event(1).domainStart=domainStart;
                    event(1).domainIncr=domainIncr;
                    event(1).startIdx=1;
                    event(1).endIdx=length(cSeq);
                    event(1).type=cSeq(1) ;
                    event(1).nSamp=event(1).endIdx-event(1).startIdx+1;
                    return
                end
                
                event(1).domainStart=domainStart;
                event(1).domainIncr=domainIncr;
                event(1).startIdx=1;
                event(1).endIdx=iChange(1);
                event(1).type=cSeq(1) ;
                event(1).nSamp=iChange(1);
                
                
                for iEvent=2:nEvent
                    event(iEvent).domainStart=domainStart;
                    event(iEvent).domainIncr=domainIncr;
                    event(iEvent).startIdx=event(iEvent-1).endIdx+1;
                    event(iEvent).endIdx=iChange(iEvent);
                    event(iEvent).type=cSeq(iChange(iEvent));
                    event(iEvent).nSamp=event(iEvent).endIdx-event(iEvent).startIdx+1;
                end
                % Last event
                event(nEvent+1).domainStart=domainStart;
                event(nEvent+1).domainIncr=domainIncr;
                event(nEvent+1).startIdx=iChange(nEvent)+1;
                event(nEvent+1).endIdx=length(cSeq);
                event(nEvent+1).type=cSeq(end);
                event(nEvent+1).nSamp=event(nEvent+1).endIdx-event(nEvent+1).startIdx+1;
                event=event';
                
                
            end
        end
        
        function event = idx2eseq ( idx, domainStart, domainIncr, eventCat )
            % idx - nx2 mstrix containing start and stop indices for each
            % event
            % domainStart
            % domainIncr
            % eventCat - category of event.
            
            [nRow,nCol]=size(idx);
            if nCol ~= 2
                error ('idx must have two columns');
            end
            if iscell(eventCat)
                eventCat=categorical(eventCat);
            elseif ~iscategorical(eventCat)
                error('eventCat must be either cstegorical or a cell array');
            end
            [nRowCat, nColCat]=size(eventCat);
            if ~ (nRowCat==1 | nRowCat ==nRow);
                error ('eventCst must be a singleton or have the same number of rows as idx');
            end
            event=eseq;
            for i=1:nRow,
                event(i).domainStart=domainStart;
                event(i).domainIncr=domainIncr;
                event(i).startIdx=idx(i,1);
                event(i).endIdx=idx(i,2);
                if nRowCat==1
                    currentCat=eventCat;
                else
                    currentCat=eventCat(i);
                end
                event(i).type=categorical(currentCat);
                event(i).nSamp=idx(i,2)-idx(i,1)+1;
                
            end
        end
    end
end

