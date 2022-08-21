function [xdot,x] = newPert(vMin,vMax,delta_v,xdotScale,minT,maxT,nSwitch,desiredpdfType,Npdf,Nbins,Nxdot,Ts)
%% This function is intended to generate the novel input perturbation signal (velocity perturbation) with a desired PDF for identification of
%% Nonlinear Hammerstein systems, to be submitted as a paper to EMBC2012 Conference

%-- Inputs: vMin and vMax: are the minimum and maximum velocities (without scaling) over which we generate the desired PDF
%           delta_v: is the velocity resolution to generate the desired PDF using histogram method (i.e., determines the width of histogram bins
%           minT and maxT: are the minimum and maximum values of pulse width (in sec) specifying random duration for the Position Perturbation Signal Pulse Widths
%           xdotScale: is the scale factor for velocity
%           nSwitch [int]: is the switching rate (in samples) between Pos and Neg velocities
%           desiredpdfType [string]: The type of velocity PDF
%           Npdf [int]: is the number of data points in the realization of the velocity PDF
%           Nbins [int]: is the number of bins to generate an estimate of the desired PDF using histogram method
%           Nxdot [int]: is the number of time samples in the generated perturbation signal
%           Ts: is the sampling time of the generated perturbation in seconds
%-- Output: xdot [Nxdot,1]: the generated velcoity perturbation
%           x [Nxdot,1]: the generated position perturbation

%% Mar. 15, 2012: This is a PDF design for perturbation generation approach
%% The PDF of desired velocity shall be biased towards large velocities and be small near zero (i.e., having a valley near zero) 
%% so that after passing through actuator low-pass dynamics, we get a relatively flat PDF for tha actual velocity. 
%++ So, in general, the PDF of desired velocity shall be V-shaped with user-controllable characteristics for the two wings of the V-shape
%++ For example, the wings can be linear, exponential, polynomial-exponential, etc.

switch desiredpdfType
    case 'V-shape'
        xdotScale = 12;  %8(pos: +/-0.015, vel: +/-1.26); %10(pos: +/-0.018, vel: +/-1.58); %12(pos: +/-0.023, vel: +/-1.9); %20; %-- Scale factor for velocity
    case 'exponential_2'
        xdotScale = 11;  %11(pos: +/-0.022, vel: +/-1.75);
end

v = vMin:delta_v:vMax;
v = v';

%% Desired PDF design 
%-- Parameters of the linear PDF. The positive half is v_pdf_p = beta + alfa_p*(v - vMin);
switch desiredpdfType
    case 'V-shape'   %++ A linear V-shape PDF
        %--- The overall idea is to have a "V-shape" distribution for ankle velocity
        %--- This means that the PDF is simply a line, the slope of which will be
        %--- determined based on range of velocity [vMin,vMax] and integral of PDF
        %--- being equal to 1 (or 0.5 for half of the symmetric PDF)
        %--- So, the parameter of PDF is slope (alfa).
        beta = 0;  %-- The intercept is beta = 0 because we want the PDF to be zero at zero velocity
        alfa_p = 1/(vMax-vMin)^2;  %-- Simple mathematical calculations wil show that the area under the positive half of the PDF is simply A = (alfa_p/2) * (vMax-vMin)^2. And the area of half PDF shall be 0.5; i.e. A = 0.5
        alfa_n = -alfa_p;
        %-- Positive half of the PDF (i.e., PDF for positive velocities)
        v_p = v;
        v_pdf_p = beta + alfa_p*v_p;
        %-- Negative half of the PDF
        v_n = -v;
        v_pdf_n = beta + alfa_n*(v_n);
        v_pdf = [v_pdf_n; v_pdf_p];
    
    case 'exponential_1'  %++  An exponential PDF
        %-- This PDF for positive half (i.e., v>=0) is PDF_p(v) = -1 + exp(alfa*v), alfa>0
        %-- If we solve for the area under PDF constraint with vMin = 0, we get the following nonlinear algebraic equation
        %-- A_p = -vMax + (exp(alfa*vMax)-1)/alfa = 0.5;
        %-- However, the above equation cannot be solved analytically and requires lambertw function 
        %-- I used MATLAB's "solve" function but didn't work! 
        %-- So, i find alfa for each vMax by plotting f(alfa) = -vMax +(exp(alfa*vMax)-1)/alfa - 0.5 and finding the point where f(alfa)~0
        alfaRange = 0:0.001:1;
        f_alfa = -vMax + (exp(alfaRange*vMax)-1)./alfaRange - 0.5;

        % figure;
        % plot(alfaRange,f_alfa);

        %++ fuinding the zero crossingf f_alfa
        j1 = find(f_alfa>=0,1,'first');
        j2 = find(f_alfa<0,1,'last');
        alfa = alfaRange(j1);
        
        %-- Positive half of the PDF (i.e., PDF for positive velocities)
        alfa_p = alfa;
        v_p = v;
        v_pdf_p = -1 + exp(alfa_p*v_p);
        %-- Negative half of the PDF
        alfa_n = -alfa_p;
        v_n = -v;
        v_pdf_n = -1 + exp(alfa_n*v_n);
        v_pdf = [v_pdf_n; v_pdf_p];
        
    case 'exponential_2'  %-- This is a PDF of the form: PDF_p(v) = alfa*(v.^3*exp(v)) ,,, In general, the int(x^n*exp(c*x)) is "(1/c)*x^n*exp(c*x) - (n/c)*int(x^(n-1)*exp(c*x))"
                          %-- For c = 1 and n= 3, we get: int(x^3*exp(x)) = exp(x)*(x^3 - 3*x^2 + 6*x - 6) 
                          %-- Doing the math for the area under positive half of the PDF, we get the following formula for alfa
        alfa = 0.5 / (6 + exp(vMax)*(vMax^3 - 3*vMax^2 + 6*vMax - 6)); 
        %-- Positive half of the PDF (i.e., PDF for positive velocities)
        v_p = v;
        v_pdf_p = alfa * (v_p.^3.*exp(v_p));
        %-- Negative half of the PDF
        v_n = -v;
        v_pdf_n = alfa * (-v_n.^3.*exp(-v_n));
        v_pdf = [v_pdf_n; v_pdf_p];        
                
    otherwise
        disp('This case is not yet designed');
        return;
end

% figure;
% subplot(2,1,1)
% plot(v_p,v_pdf_p);
% title('Positive half of the PDF');
% subplot(2,1,2)
% plot(v_n,v_pdf_n);
% title('Negative half of the PDF');

% figure;
% plot([v_n;v_p],v_pdf,'*')


%% To generate a realization of a random sequence with the desired PDF
%-- I would use a binning approach
% Ts = 0.001;   %--- Sampling time
time = 0:Ts:(Npdf-1)*Ts;
vMax_p = vMax;  %-- rad/sec
vMin_n = -vMax_p;  %-- rad/sec
deltaV = (vMax_p - vMin_n)/Nbins;
vBins = vMin_n:deltaV:vMax_p;
nPerBin = zeros(Nbins,1);

vPert = zeros(Npdf,1);

%+++ For the "V-shape" PDF
i1 = 1;
for i = 1:Nbins
    if vBins(i) < 0 
        switch desiredpdfType
            case 'V-shape'
                h1 = beta + alfa_n*vBins(i);
                h2 = beta + alfa_n*vBins(i+1);
            case 'exponential_2'
                h1 = alfa * (-vBins(i).^3.*exp(-vBins(i)));
                h2 = alfa * (-vBins(i+1).^3.*exp(-vBins(i+1)));
        end
    else
        switch desiredpdfType
            case 'V-shape'
                h1 = beta + alfa_p*vBins(i);
                h2 = beta + alfa_p*vBins(i+1);
            case 'exponential_2'
                h1 = alfa * (vBins(i).^3.*exp(vBins(i)));
                h2 = alfa * (vBins(i+1).^3.*exp(vBins(i+1)));
        end
    end
    Area = 0.5*(h1+h2)*deltaV;
    nPerBin(i) = fix(Area*Npdf);
    i2 = i1 + nPerBin(i) - 1;    
    vPert(i1:i2) = vBins(i) + (vBins(i+1)-vBins(i))*rand(nPerBin(i),1);
    i1 = i2 + 1;
end

%% Randomizing the realization sequence with the desired PDF
randInd = randperm(Npdf);
vPertRand = zeros(Npdf,1);
for j = 1:Npdf
    vPertRand(j) = vPert(randInd(j));
end

% figure;
% plot(time,vPertRand)

% vPertRand_pdf = pdf(vPertRand,'NBins',100,'pdfType','Density');
% figure;
% plot(vPertRand_pdf)

%% Generating the velocity signal samples
time = 0:Ts:(Nxdot-1)*Ts; time = time';
cnt1 = 1;
cnt2 = 1;
xdot = -100*ones(Nxdot,1);
xdotImpulse = zeros(fix(Nxdot/minT/Ts),1);
% xdotImpulseNeg = -xdotImpulse; 
while cnt1 < Nxdot
    %-- Every "nSwitch" samples of the xdotImpulse, we switch to negative of velocity Impulses
    n = floor((cnt2-1)/nSwitch);
    if rem(n,2)==0    %--- Checking if quotient of cnt2 & nSwitch is ODD or EVEN
        %-- Picking up a random velocity impulse from the relaization of the V-Shape distribution
        randIndices = randperm(Npdf);
        randInd = randIndices(1);
        xdotImpulse(cnt2) = xdotScale*vPertRand(randInd);
        %-- Generating a random duration for zero velocity 
        T = fix(minT + (maxT - minT)*rand);
        xdotZero = zeros(T,1);
        xdotPulse = [xdotImpulse(cnt2);xdotZero];
        cnt2 = cnt2 + 1;
        cnt3 = 1;
    else   %--- Now, we should use the negative of xdotImpulse in a randomized fashion
        if cnt3 == 1 %rem(cnt2,nSwitch) == 1
            xdotImpulseNeg = -xdotImpulse(cnt2-nSwitch:cnt2-1,1);
            randIndices = randperm(nSwitch);
            xdotImpulseNegRand = xdotImpulseNeg(randIndices);
        end
        %-- Generating a random duration for zero velocity 
        T = fix(minT + (maxT - minT)*rand);
        xdotZero = zeros(T,1);
        xdotPulse = [xdotImpulseNegRand(cnt3);xdotZero];
        cnt2 = cnt2 + 1;
        cnt3 = cnt3 + 1;
    end
    %-- Assigning to xdot
    i1 = cnt1;
    i2 = i1 + length(xdotPulse) - 1;  %-- length(xdotPulse) = T+1
    xdot(i1:i2) = xdotPulse;
    cnt1 = i2 + 1;
end

xdot = xdot(1:Nxdot,1);
x = cumsum(xdot)*Ts;

% NFFT = 5000;

% xdot_pdf = pdf(xdot,'NBins',100,'pdfType','Density');
% xdotData = nldat(xdot,'domainIncr',Ts);
% xdot_spect = spect(xdotData,'NFFT',NFFT);

% figure;
% subplot(3,1,1)
% plot(time,xdot); title('Velocity Perturbation');
% xlim([40,50]);
% xlabel('time (sec)')
% ylabel('velocity (rad/s)')
% subplot(3,1,2)
% plot(xdot_pdf)
% ylim([0,0.002])
% subplot(3,1,3)
% plot(xdot_spect)

% x_pdf = pdf(x,'NBins',100,'pdfType','Density');
% xData = nldat(x,'domainIncr',Ts);
% x_spect = spect(xData-mean(xData),'NFFT',NFFT);

% figure;
% subplot(3,1,1)
% plot(time,x); title('Position Perturbation');
% xlim([40,50]);
% xlabel('time(sec)')
% ylabel('position (rad)')
% subplot(3,1,2)
% plot(x_pdf)
% % ylim([0,0.02])
% subplot(3,1,3)
% plot(x_spect)









