function [ifil] = invert_irf(fil,dt,invf_len);
%  Inverts an impulse response function in the frequency domain
%  Both input and output IRFs are assumed to be
%  two-sided and symmetric.
%  
%  ifil=invert_irf(fil,dt,invf_len) 
%

fftlen=(fix(log(invf_len)/log(2)));  % Choose FFT length to be nearest
if rem(2^fftlen,invf_len)~= 0
    fftlen=fftlen+1;                  % power of 2 > invf_len  
end
fftlen=2^fftlen

padr=0;
fil=fil(:);

%  Add both leading and trailing zeroes to IRF to meet necessary FFT length:
if fftlen > length(fil)
   pad=fftlen-length(fil)
   pad2=pad/2;
   padr=rem(pad,2)
   fil=[zeros(pad2,1);fil;zeros(pad-pad2,1)];
end

%  Take FFT, take its reciprocal, then invert via inverse FFT:
   filfft=fft(fil,fftlen)*dt;
   ifilfft=1.0 ./filfft;
   ifil=real(ifft(ifilfft,fftlen)/dt);

%  FFT length will always be even, while IRF length will be odd;  extra
%  value to delete is the first sample:
   ifil=ifil(2:length(ifil)); 

return



