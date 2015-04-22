function bl = phk_ident(bl,iodata);
% implements Paulin-Hunter-Korenberg interation for Wiener cascade models.
% 
% See M.G. Paulin, A method for constructing data-based models of spiking 
%     neurons using a dynamic linear-static nonlinear cascade.
%     Biological Cybernetics, 69:67-76, 1993.
%
%     M.J. Korenberg and I.W. Hunter,  Two methods for identifying Wiener 
%     cascades having noninvertible static nonlinearities.
%     Annals of Biomedical Engineering, 27(6):793-804, 1999.

% Copyright 2003, Robert E Kearney and David T Westwick
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
ErrGain = get(bl,'updateGain');

x_est = nlsim(h_est,u);
y_est = nlsim(m_est,x_est);
err = z - y_est;
old_vaf = double(vaf(z,y_est));

count = 1;
while count < MaxIts
  x_est2 = x_est + ErrGain * err;
  uxe2 = cat(2,u,x_est2);
  h_new = nlident(h_est,uxe2);
  x_est = nlsim(h_new,u);
  xez = cat(2,x_est,z);
  m_new = nlident(m_est,xez);
  y_est = nlsim(m_new,x_est);
  err = z - y_est;
  new_vaf = double(vaf(z,y_est));
  if new_vaf - old_vaf > ErrTol
    h_est = h_new;
    m_est = m_new;
    old_vaf = new_vaf;
    count = count + 1;
  else
    count = MaxIts + 1;
    if new_vaf > old_vaf
      h_est = h_new;
      m_est = m_new;
    end
  end
end

subsys = {h_est,m_est};
set(bl,'elements',subsys,'comment',...
    'LN model identified using Paulin-Hunter-Korenberg iteration');




