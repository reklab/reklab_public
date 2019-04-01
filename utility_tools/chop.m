function c = chop(x, t)
%CHOP    Round matrix elements.
%        CHOP(X, t) is the matrix obtained by rounding the elements of X
%        to t significant decimal places.

scale=10^t;
c = fix(x*scale)/scale;