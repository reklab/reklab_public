 
function [B,D,x0,K,R,Phi] =  cestbd(u,y,Ts,A,C,model,Rold); 
  
% cestbd   Estimates the matrices B and D of the state space model 
% 
%                  . 
%                  x(t) = Ax(t) + Bu(t) + w(t) 
%                  y(t) = Cx(t) + Du(t) + v(t) 
%           using the knowledge of the pair A, C. This function enables 
%           to concatenate different input-output data batches, through 
%           the matrix R and Rold). 
% 
%           B and D and x0 are calculated by solving a linear least squares 
%           problem. 
% 
% 
% Syntax: 
%           [B,D,x0,R,Phi]=cestbd(u,y,Ts,A,C,[fB fD fx],Rold); 
%           [B,D]=cestbd(u,y,Ts,A,C); 
% 
% Input: 
%   u, y    The input respectively output data of the system to be 
%           identified. 
%   Ts      Sampling period of the measured data. 
%   A, C    The estimated system matrices A and C. 
%   model   Three element flag vector [fB fD fx] indicating whether 
%           B, D and x0 should be estimated. The default value is [1 1 0]. 
%           The matrix B or D  can be assumed zero by setting 
%           fB or fD to zero. The calculation of x0 can be ommitted by 
%           setting fx to zero. x0 will not be assumed zero in 
%           this case It's influence will still be taken into account for the 
%           computation of B and D. However, not calculating x0 makes 
%           the calculation more robust. 
%   Rold    R matrix obtained from previous data batch. This variable 
%           can be used to process data in batches, or to combine data 
%    from different experiments. 
% 
% Output: 
%   B,D     The estimated B, D pair. 
%   x0      Estimated initial condition. 
%   R       Compressed data matrix, storing information on the calcuation of 
%           the matrices B and D in following cestbd. Used when analyzing 
%           multiple input-output data sequences. 
% 
% see also: cestx 
 
%  --- This file is generated from the MWEB source cestbd.web --- 
% 
% Bert Haverkamp, August 1999 
% copyright (c) 1999 Bert Haverkamp, Michel Verhaegen 
 
 
 
if nargin==0 
  help cestbd 
  return 
end 
if nargin<6 
  model = []; 
end 
if isempty(model) 
  model = [1 1 0 0]; 
elseif length(model)==3 
  model = [model,0]; 
elseif length(model)~=4 
  error('Variable ''model'' should have four elements') 
end 
 
if nargin<7 
  Rold = []; 
end 
 
fB = model(1); 
fD = model(2); 
fx = model(3); 
fK = model(4); 
I  =  eye(length(A)); 
[Lt,Ut,Pt]  =  lu(I - A * Ts/2); 
Ad  =  (((I + A * Ts/2)/Ut)/Lt) * Pt;  % (I + a*T/2)/(I - a*T/2) 
Cd  =  Ts * ((C/Ut)/Lt) * Pt;          % T*c/(I - a*T/2) 
[Bd,Dd,x0d,Kd,R] = destbd(u,y,Ad,Cd,model,Rold); 
[Lt,Ut,Pt]  =  lu(I + Ad); 
if fB 
  B  =  2 * (Ut \(Lt \(Pt * Bd))); 
else 
  B = []; 
end 
if fD 
  if fB 
    D  =  Dd - Cd * (B/2); 
  else 
    D = Dd; 
  end 
else 
  D = []; 
end 
if fx 
  x0 = x0d * Ts; 
else 
  x0 = []; 
end 
if fK 
  K  =  2 * (Ut \(Lt \(Pt * Kd))); 
else 
  K = []; 
end 
 

