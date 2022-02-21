%% This function generates Laguerre orthonormal basis functions
%++ Author: Ehsan Sobhani (10 April 2014)
%++ This is based on Maremaleris book OR formula (11) of his paper titled:
%++ "Identification of Nonlinear Biological Systems Using Laguerre Expansions of Kernels", Annals of Biomed. Eng., vol. 21, pp. 573-589, 1993.

function b = laguerre(irf_len,max_order,alfa)

%++ The inputs are:
    % 1) irf_len (integer number of samples)
    % 2) max_order (integer). The maximum order of Laguerre expansion.
    % 3) alfa (0<alfa<1). Laguerre parameter. 

%++ The outputs are:
    % 1) b (real with size [irf_len,max_order+1])

%=== Initialization
b = zeros(irf_len,max_order+1);
    
for t = 1:irf_len
    for j = 0:max_order
        gain = alfa^((t-j)/2) * (1-alfa)^(0.5);
        summation = 0;
        for k = 0:j
            argument = (-1)^k * combination(t,k) * combination(j,k) * alfa^(j-k) * (1-alfa)^k;
            summation = summation + argument;
        end
        b(t,j+1) = gain * summation;
    end
end    

end

