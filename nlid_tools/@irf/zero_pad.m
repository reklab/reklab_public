function IRF_Padded = zero_pad(IRF_In);
% new_IRF_obj = zero_pad(IRF_obj)
%
% IRF_obj: tv_IRF object
%
% Pad a TV IRF matrix(and corresponding bounds matrix if not empty) 
% with M2 rows of zeros at the top and -M1 rows of zeros at the bottom.
% so that it can be used efficiently for simulation
%

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

P=IRF_In.parameterSet;
assign(P);
ds=IRF_In.domainstart;
dt=IRF_In.domainincr;
%
% Determine padding needed
%
hlength = size(IRF_In,2);
if nSides==1,
    M2=hlength;
    M1=0;
else
    M2=(hlength-1)/2;
    M1=-M2;
end

% Pad the IRFs matrix.  
hlength = size(IRF_In,2);
IRF_Data= IRF_In.dataSet;
IRF_Data=squeeze(IRF_Data);
IRF_Data = [zeros(M2,hlength); IRF_Data; zeros(-M1,hlength)];
ds(1)=0;

% % If the bounds matrix is not empty, pad the bounds matrix too.
% if ~(isempty(bounds))
%    bounds = [zeros(M2,hlength); bounds; zeros(-M1,hlength)];
% end

% Create the new IRF object.
IRF_Padded = IRF_In;
set(IRF_Padded,'dataSet',IRF_Data,'domainStart',ds);


% end @tv_IRF\zero_pad.m
