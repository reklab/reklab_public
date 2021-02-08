function event =cseq2eseq(cSeq, caseNum)
%catSignal2Events - convert a categorical sequence  to event sequence
% Input: cSig - categorical signal
% output 
% event - structure with fields
% event.type 
% event.start
% event.end,
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

event(1).start=1;
event(1).end=iChange(1);
event(1).type=cSeq(1) ;
event(1).length=event(1).end-event(1).start+1; 
event(1).caseNum=caseNum; 

for iEvent=2:nEvent  
    event(iEvent).start=event(iEvent-1).end+1;
    event(iEvent).end=iChange(iEvent);
    event(iEvent).type=cSeq(iChange(iEvent));
    event(iEvent).length=event(iEvent).end-event(iEvent).start+1;
    event(iEvent).caseNum=caseNum;

end

event(nEvent+1).start=iChange(nEvent)+1;
event(nEvent+1).end=length(cSeq);
event(nEvent+1).type=cSeq(end);
event(nEvent+1).length=event(nEvent+1).end-event(nEvent+1).start+1; 
event(nEvent+1).caseNum=caseNum;
event=event';
end









