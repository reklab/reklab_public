function compute_pinv(control);
%  work-horse routine for pseudo-inverse based IRF computation.
%  syntax:   compute_pinv(control);
%
%  control is either 'init', in which case initial correlation
%  computations are performed, or "irf", in which case an irf and confidence
%  bound are computed (of order dtw_pinv_order)
%
%
%  Operations are performed on a global workspace.  
%  The following globals are assumed to have been created prior to
%  compute_pinv being called:
%
%  for option 'init'
%
%  dtw_pinv_input dtw_pinv_output dtw_pinv_numlags dtw_pinv_numsides
%
%  for option 'irf', 'init' must have been called with the current data, and
%  the following must also be defined:
%
%  dtw_pinv_order
%
%  the golbal variables dtw_pinv_irf and dtw_pinv_bound are used to return
%  the irf and bound, if option "irf" is selected.

% Copyright 1994-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

global dtw_pinv_input dtw_pinv_output dtw_pinv_numlags dtw_pinv_numsides
global dtw_pinv_order dtw_pinv_level

global dtw_auto dtw_yy dtw_cross dtw_U dtw_V dtw_Sinv dtw_mdls
global dtw_out_var dtw_rand_lev dtw_bias_lev dtw_phi_uu2

global dtw_pinv_irf dtw_pinv_bound

N = length(dtw_pinv_input);
numlags = dtw_pinv_numlags;
numsides = dtw_pinv_numsides;

if strcmp(control,'init')
  dtw_pinv_input = detrend (dtw_pinv_input);
  dtw_pinv_output = detrend (dtw_pinv_output);
  [dtw_auto,dtw_yy,dtw_cross] = ...
      scorr(dtw_pinv_input,dtw_pinv_output,dtw_pinv_numlags); 
  if numsides == 1
    dtw_auto = dtw_auto(numlags:2*numlags-1);
    dtw_cross = dtw_cross(numlags:2*numlags-1);
  else
    dtw_auto = dtw_auto(numlags:2*numlags-1);
    i1=round((numlags+1)/2);
    i2=round((3*numlags-1)/2);
    dtw_cross = dtw_cross(i1:i2);
  end
  dtw_auto = dtw_auto(:);
  dtw_cross = dtw_cross(:);

  [dtw_U,S,dtw_V] = svd(toeplitz(dtw_auto));
  S = diag(S);
  dtw_Sinv = diag(1./S);
  
  %  sort the "modes" in order of their output contribution, rather than in
  %  order of input power (default provided by the svd)
  
  out_var = dtw_Sinv*((dtw_U'*dtw_cross).^2);
  [coeff,indeces] = sort(out_var);
  indeces = flipud(indeces);
  coeff = flipud(coeff);
  dtw_V = dtw_V(:,indeces);
  dtw_U = dtw_U(:,indeces);
  S = S(indeces);
  dtw_Sinv = diag(1./S);

  % fast calculation of mdls   
  dtw_mdls = N + [1:numlags]'*log(N);
  dtw_mdls = dtw_mdls.*(dtw_yy(numlags) - cumsum(coeff));

  %  Since the sorting based on output power was incorporated, the
  %  fast mdl computation was added "inline", since some of the 
  %  computations had already been done -- the whole calculations
  %  can be performed by the following:
  %dtw_mdls = fast_mdl1(dtw_yy(numlags),dtw_Sinv,dtw_U,dtw_cross,N);

  
  %  common calculations required in all bias error estimates
  h0 = toep(dtw_auto,dtw_cross);
  resid = dtw_pinv_output -  filter_ts(h0,dtw_pinv_input,numsides);
  coeff = dtw_V'*h0;
  dtw_out_var = cumsum((coeff.^2).*S);
  dtw_out_var = dtw_out_var/dtw_out_var(numlags);
  %  Determine the confidence level for the bias component.
  %  Tchebyshev polynomial from bias confidence level to 
  %  actual bound
  poly_clev = [1; 99; 0.59689110450194; 0.35661622346423;...
                      0.05302874116126; 0.07068789485648;...
                      0.02058722980738; 0.02726886587237;...  
                      0.01013416651667; 0.01135442050998];
  dtw_bias_lev = tchebval(poly_clev,dtw_pinv_level); 


  % common calculations required in all random error estimates
  dtw_phi_uu2 = [flipud(dtw_auto); dtw_auto(2:numlags)];
  %  generate a Gaussian cdf  note that we are only generating the 
  %  positive half, and then normalizing wrt the sum over that half,
  %  so the cdf correponds to the probability of being between the limits 
  dom = [0:0.01:5];
  gauss_pdf = exp(-(dom.^2)/2)/sqrt(2*pi);
  gauss_cdf = cumsum(gauss_pdf)/sum(gauss_pdf);
  dtw_rand_lev = dom(min(find(gauss_cdf>(dtw_pinv_level/100))));



elseif strcmp(control,'irf')
  order = dtw_pinv_order;
  PHIuui = dtw_V(:,1:order)*dtw_Sinv(1:order,1:order)*dtw_U(:,1:order)';
  dtw_pinv_irf = PHIuui*dtw_cross;
 
  testout = filter_ts(dtw_pinv_irf,dtw_pinv_input,numsides);
  resid = dtw_pinv_output - testout;
  bias_norm = norm(dtw_pinv_irf)*sqrt(1-dtw_out_var(dtw_pinv_order));
  phivv = phixy(resid,numlags);
  phi_vv2 = [flipud(phivv); phivv(2:numlags); zeros(2*numlags-1,1)];
  test = filter(dtw_phi_uu2,1,phi_vv2);
  test = test(2*numlags-1:3*numlags-2);
  PHIuv = toeplitz(test)/N;
  hvar = abs(diag(PHIuui*PHIuv*PHIuui));
  dtw_pinv_bound = sqrt((dtw_rand_lev^2)* hvar + (dtw_bias_lev*bias_norm)^2);
end  


