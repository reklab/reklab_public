function cseq = eseq2cseq (eseq)
% eseq2cseq - converts an event sequent to to a categorical sequence

cseq=categorical;
for i=1:length(eseq)
  cseq(eseq(i).start:eseq(i).end,1)=eseq(i).type; 
end


