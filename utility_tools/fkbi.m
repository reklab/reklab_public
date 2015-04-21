function [y,dyda]=fkbi(x,a)
% compute second-order IRF for K,B,I
%
% Function to fit K, B, and I of a 2nd order compliance
% to an impulse response

% Usage: [y,dyda]=fkbi(x,a)
%                x    ==  lag (s)
%          where y       ==    compliance IR at time x
%                a=[k b i ]
%                   k  ==   elastic stiffness (Nm/rad)
%                   b  ==   viscous stiffness (Nm/rad/s)
%                   i  ==   inertia (Nm/rad/s^2)
%                dIda    ==    matrix of partials
% 21 Feb 95 REK make dyda real.

x=x(:);  
k=a(1);
b=a(2);
i=a(3);
t=x;

	dyda=zeros(length(x),length(a));

	c1=t*sqrt((k/i)-((b^2)/(4.*i*i)));
	c2=exp((b*t)/(2.*i));
	c3=(k/i)-((b^2)/(4.*i*i));
	c4=(0.5*b*b)/(i^3)-k/(i*i);

	dyda(:,1)=(0.5*t.*cos(c1))./(c2*i*i*c3) - ...
       ((0.5*sin(c1))./(c2*i*i*(c3^1.5)));

        dyda(:,2)=(-0.25*b*t.*cos(c1))./(c2*i*i*i*c3) + ...
       (0.25*b*sin(c1))./(c2*i*i*i*(c3^1.5)) - ...
       (0.5*t.*sin(c1))./(c2*i*i*sqrt(c3));

        dyda(:,3)=(0.5*t.*cos(c1)*c4)./(c2*i*c3) - ...
       (0.5*c4*sin(c1))./(c2*i*(c3^1.5)) - ...
       (sin(c1))./(c2*i*i*sqrt(c3)) + ...
       (0.5*b*t.*sin(c1))./(c2*i*i*i*sqrt(c3));


        e=(1.)/(i*sqrt(c3));
	f=(1.)./c2;
	g=sin(c1);
	
	y=e.*f.*g;




% The partial derivative equations generate "intermediate" complex results
% which combine to give a real results, but due to round off error, the
% derivatives are often of the form: real + 0.000000i.  The identification
% function (nlls) works fine even when the derivatives have this form,
% (I have verified this)  but since I know the result is always real, I will
% due away will the extraneous imaginary part, also, I reshape dcda so that
% it works properly with nlls:
dyda = real(dyda);


% end of file fkbi.m
