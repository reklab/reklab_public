function nlmtst(w)
% nlmtst for waveforms
t=0:.001;1;
w=waveform('waveformType','pulse','delay',.1, 'width',.5);
y=nlsim(w,t);
plot(y)
