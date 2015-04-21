 function [A,B,C,D,K] = dmoesp(u,y,maxorder,order,method) 
  
% dmoesp  High level function for the discrete time moesp that gives 
%         the user a simple interface to the lower level moesp 
%         functions such as  dordpo and destac.  The function 
%         calculates the system matrices A, B, C and D from the 
%         given input and output data sequences u and y. 
%         Model structure: 
%            x(k+1) = Ax(k) + Bu(k) 
%            y(k)   = Cx(k) + Du(k) 
%         For information about the possible disturbance signals see 
%         the help pages for the preprocessor {\dordxx} functions. 
% 
%         The user only needs to specify the input and output data. 
%         If the order is not given, a order selection dialog is 
%         shown where the user can select the order by examining 
%         the singular value plot. 
% 
% Syntax: 
%           [A,B,C,D]=dmoesp(u,y,maxorder) 
%           [A,B,C,D,K]=dmoesp(u,y,maxorder,order,method) 
% 
% Input: 
% u,y       Input and Output data sequence. 
% maxorder  Maximum expected order of the system, this value is 
%           used as the Hankel dimension parameter in the underlying 
%           identification routines. 
% order     Order of the system to be estimated, if not specified a 
%           dialog will be presented which allows the user to choose the 
%           order. 
% method    Method used for identification. 
%           possible options: 'om'  = ordinary moesp 
%                             'pi'  = past input moesp 
%                             'po'  = past output moesp 
%           The default value is past output moesp. 
% 
% Output: 
% A,B,C,D   State space matrices describing the estimated model. 
% K         Estimated Kalman gain, only available for past output moesp. 
% 
% See also: destac dordpo, etc 
 
%  --- This file is generated from the MWEB source dmoesp.web --- 
% 
% Bert Haverkamp 12/11/99 
% Copyright (c) 1999 Bert Haverkamp 
 
 
 
 if nargin==0 
  help dmoesp 
  return 
end 
argmin = 3; 
 if nargin<argmin 
  error('There are not enough arguments.') 
end 
 
 
 
if nargin <4 
  order = []; 
end 
 
if isempty(order) 
  ordersel_flag = 1; 
else 
  ordersel_flag = 0; 
end 
 
if nargin <5 
  method = 'po'; 
end 
 
if length(method)>2 
  method = method(1:2); 
end 
 
if (method=='po')& (nargout>4) 
  fK = 1; 
else 
  fK = 0; 
end 
 
 
 
 
 
  
 
if method=='om' 
  [S,R]  =  dordom(u,y,maxorder); 
  if ordersel_flag 
    order = orderselect(S); 
    if isempty(order) 
      A = [];B = [];C = [];D = []; 
      return; 
    end 
  end 
elseif method=='pi' 
  [S,R]  =  dordpi(u,y,maxorder); 
  if ordersel_flag 
  order = orderselect(S); 
    if isempty(order) 
      A = [];B = [];C = [];D = []; 
      return 
    end 
  end 
elseif method=='po' 
  [S,R]  =  dordpo(u,y,maxorder); 
  if ordersel_flag 
    order = orderselect(S); 
    if isempty(order) 
      A = [];B = [];C = [];D = [];K = []; 
      return 
    end 
  end 
else 
  error('unknown method') 
end 
[A,C]  =  destac(R,order); 
[B,D]  =  destbd(u,y,A,C); 
if fK 
  K = destk(A,C,R); 
end 
  
 
 
 

