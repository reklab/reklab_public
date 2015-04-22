function [y,x]=tlsim(A,B,C,D,u,Ts,x0)
% tlsim  Simulation of the time response of LTI continuous-time 
%        systems to arbitrary inputs. 
%        The input and simulated output are represented by 
%        sampled signals. Unlike lsim, where the sampling is assumed 
%        to be zero order or first order hold, the sampling in tlsim is
%        assumed to be trapezoidal or Tustin's (hence the t in tlsim).  
%
% syntax
%  [y,x]=tlsim(A,B,C,D,u,Ts,x0)
%
% Input:
%  A,B,C,D  State space matrices of the LTI system.  
%  u        Input signal
%  Ts       Sampling time.
%  x0       Initial state of the system.
%  
%  Output:
%  y        Output of the LTI system.
%  x        State of the LTI system. 
%
% See also: lsim, dlsim  

if nargin==0
  help tlsim
  return
end
if length(Ts)>1
  Ts=Ts(2)-Ts(1);
end
if nargin<7
  x0=zeros(length(A),1);
end
I  =  eye(size(A)); 
[L,U,P]  =  lu(I - A * Ts/2);
Ad  =  (((I + A * Ts/2)/U)/L) * P;   % (I + a*T/2)/(I - a*T/2) 
Cd  =  Ts * ((C/U)/L) * P;  % T*c/(I - a*T/2)
Bd = U\(L\(P*B));              % (I - a*T/2)\b
Dd = Cd*B/2 + D;        % (T/2)*c*(I - a*T/2)\b
x0d=x0/Ts;              
[y,x]=dlsim(Ad,Bd,Cd,Dd,u,x0d);

%sysc=ss(A,B,C,D);
%sysd=c2d(sysc,Ts,'tustin')
