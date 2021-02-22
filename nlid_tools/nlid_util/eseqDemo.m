function eseqDemo 
% eseqDemo - domnonsrate properties of the event sequence class

e1= eseq;
for i=1:5;
    e1(i,1).startIdx=1+i*20;
    e1(i,1).endIdx=e1(i).startIdx+10;
    e1(i,1).type='1';
    
    e2(i,1).startIdx=1+(i-1)*20;
    e2(i,1).endIdx=e2(i).startIdx+10;
    e2(i,1).type='1';
    
    
end

% convert to cseq
cseq(e1);

% plot cseq 
figure(1);
demoPlot(e1,e2); 
 figure(2);

e2(1).startIdx=5;
demoPlot(e1,e2); 

%%
figure(3);
clf

e2(3).startIdx=e1(3).startIdx+5;
e2(3).endIdx=e2(3).startIdx+10;
demoPlot(e1,e2);
disp('tes');




end

function demoPlot (e1,e2)
e3=intersect(e1,e2);
subplot (3,1,1); plot (e1);
subplot (3,1,2); plot (e2); 
subplot (3,1,3); plot (e3); 
end





