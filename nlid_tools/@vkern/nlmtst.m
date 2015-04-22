function vk=nlmtst(vk)
%
% test of vkern objects%
% vkernal - fast orthogoal
z=rand(10,10);
vk=vkern(z,'kernOrder',2);
plot(vk);

x=nldat(rand(1000,1));
y=nlsim(vk,x);
end