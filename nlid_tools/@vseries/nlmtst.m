function vs=nlmtst(i)
%
% test of vseries identification
%
z=nlid_sim ('LN2');
% vseries - fast orthogoal
vs=vseries(z,'nLags',20);
figure(1);
plot(vs);
zp=nlsim(vs,z(:,1));
figure(2)
nlid_resid ( vs, z);

% vkern/nlmtst
