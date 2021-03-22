function [ApEn] = ApEn(x)
% Approximate Entropy of signals
% Input: 
%       - x: nldat object a signal and sampling period
%       information
% Outputs:
%       - ApEn: Approximate Entropy of the HR signal
% ----------------------------------------------------------------
% Johann Vargas-Calixto 03-10-2021

% Initialize variables
x=double(x);

% Removes useless samples
x(isnan(x))=[];

% Estimate ApEn
ApEn=approximateEntropy(x);

end

