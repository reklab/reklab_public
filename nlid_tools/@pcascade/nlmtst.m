function pc=nlmtst(i)
%
% test of nlbl identification
%

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 
clear all
x=nldat(randn(10000,1),'domainIncr',.001);
y= pcModel(x);
z=cat(2,x,y);
pc=pcascade;
pc=nlident(pc,z,'idMethod','eig');
%[r,v]=nlid_resid(pc,z); plot(pc)
%pc=nlident(pc,z,'idMethod','slice');[r,v]=nlid_resid(pc,z); plot(pc)
%pc=nlident(pc,z,'idMethod','sls');
%pc=nlident(pc,z,'idMethod','lm');[r,v]=nlid_resid(pc,z); plot(pc)
%pc=nlident(pc,z,'idMethod','gen_eig');





% @pcascade/nlmtst
