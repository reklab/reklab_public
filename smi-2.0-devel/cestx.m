 function x0 = cestx(u,y,Ts,A,B,C,D); 
  
% cestx     Estimate the initial state, given the estimated system 
%           matrices, and a set of input/output data. 
% 
% Syntax: 
%           x0=cestx(u,y,Ts,A,B,C,D); 
% 
% Input: 
%   u, y    The input respectively output data of the system to be 
%           identified. 
%   Ts      Sampling period of the measured data. 
%   A,B,C,D System matrices. 
 
% 
% Output: 
%   x0      Estimated initial state. 
% 
% See also: cestbd 
 
%  --- This file is generated from the MWEB source cestbd.web --- 
% 
% Bert Haverkamp, April 1996 
% copyright (c) 1996 B.R.J. Haverkamp 
 
 
 
 
if nargin==0 
  help cestx 
  return 
end 
I  =  eye(length(A)); 
[Lt,Ut,Pt]  =  lu(I - A * Ts/2); 
Ad  =  (((I + A * Ts/2)/Ut)/Lt) * Pt;  % (I + a*T/2)/(I - a*T/2) 
Bd  =  Ut \(Lt \(Pt * B));             % (I - a*T/2)\b 
Cd  =  Ts * ((C/Ut)/Lt) * Pt;          % T*c/(I - a*T/2) 
Dd  =  Cd * B/2 + D;                 % (T/2)*c/(I - a*T/2)*b 
x0d = destx(u,y,Ad,Bd,Cd,Dd); 
x0 = x0d * Ts; 
 

