function bl = hk_ident(bl,iodata);
% implements Hunter-Korenberg interation for Wiener cascade models.
% 
% See: I.W. Hunter and M.J. Korenberg, The identification of nonlinear 
%      biological systems: Wiener and Hammerstein cascade models.
%      Biological Cybernetics, 55:135-144, 1986.

% Copyright 1991-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


subsys = get(bl,'elements');
h_est = subsys{1};
m_est = subsys{2};

h_mode = get(h_est,'irfPseudoInvMode');
m_mode = get(m_est,'polyOrderSelectMode');


u = iodata(:,1);
z = iodata(:,2);

MaxIts = get(bl,'nMaxIts');
ErrTol = get(bl,'iterationTolerance');


x_est = nlsim(h_est,u);
zxe = cat(2,z,x_est);
minv = nlident(m_est,zxe);
x_est2 = nlsim(minv,z);
old_vaf = double(vaf(x_est,x_est2));

count = 1;
while count <= MaxIts
  uxest2 = cat(2,u, x_est2);
  h_new = nlident(h_est, uxest2);
  x_est = nlsim(h_new,u);
  zxe = cat(2,z,x_est);
  minv = nlident(m_est,zxe);
  x_est2 = nlsim(minv,z);
  new_vaf = double(vaf(x_est,x_est2));
  if new_vaf - old_vaf > ErrTol
    % improvement is significant, so continue
    h_est = h_new;
    old_vaf = new_vaf;
    count = count + 1;
  else
    count = MaxIts + 1;
  end
end

x_est = nlsim(h_est,u);
xez = cat(2,x_est,z);
set(m_est,'polyOrderSelectMode',m_mode);
set(h_est,'irfPseudoInvMode',h_mode);
m_est = nlident(m_est,xez);
subsys = {h_est,m_est};
set(bl,'elements',subsys,'comment',...
    'LN model identified using Hunter-Korenberg iteration');




