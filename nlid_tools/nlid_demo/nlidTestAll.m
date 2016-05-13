function [ output_args ] = nlid_test( input_args )
% nlid_test - test functionality of all nlid objects.
%   Detailed explanation goes here
dbstop if error
classList= { 'nldat' 'kern' 'cor', 'irf' 'spect' 'fresp' 'pdf' 'polynom' ...
    'nlbl' 'lnbl' 'pcascade' 'versies' 'wseries'};
for iClass = 1:length(classList),
    curClass=classList{iClass};
    disp([ 'Start Testing: ' curClass]);
    nlmtst(eval(curClass));
     disp([ 'Done Testing: ' curClass]);
end




end

