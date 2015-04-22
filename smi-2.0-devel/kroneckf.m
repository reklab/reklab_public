function [AA,BB,Q,Z,na,nb]=kroneckf(A,B,alpha,option);

% kroneckf  Calculates the Kronecker canonical form from a regular 
%           pencil A, B such that 
%                                            | Ak    0 |
%                                 Q*A*Z=AA=  |       I |
%          and 
%                                 Q*B*Z=BB=  | I     0 |
%                                            | 0     Bk|
%          where the eigenvalues of As and Bs are split according to
%          alpha and option.
% Syntax:
%          [AA,BB,Q,Z,na,nb]=kroneckf(A,B,alpha,option)
%
% Input:
%  A,B     The given matrix pencil.
%  alpha   The value used to split the eigenvalues of the pencil A and B.
%          Default value is 1.
%  option  This parameter can be either 'circle' or 'halfplane'.
%          'circle' sorts the eigenvalues in Ak to be absolutely
%          smaller than alpha,  and those in Bk to be larger than alpha
%          'halfplane' sorts the eigenvalues in Ak to have smaller real
%          part than alpha, and those in Bk to have larger real part.
%         
% Output:
%  AA,BB   The pencil in canonical form.
%  Q       The required left transformations.
%  Z       The required right transformations.
%  na      Dimension of the matrix Ak.
%  nb      Dimension of the matrix Bk.
%
% Seealso ncdestac 

%    Bert Haverkamp 4-5-95
%    copyright (c) 1995 Bert Haverkamp

%%%%%%% Argument check
if nargin==0
  help kroneckf
  return
end
if nargin==2  
  alpha=1;option='circle';
end

if strcmp(option,'circle') & alpha==0
  error('If option is ''circle'', alpha must be larger than 0')
end

if size(A)~=size(B)
  error('A and B should have the same size')
end
beta=eps*norm(A)*norm(B);
n=length(A);

%%%%%% Cowcatcher
if 0
  %abs(eig(A,B))==beta*ones(n,1)
  fprintf('Matrix pensil contains all zero eigenvalues.\n');
  fprintf('I can''t calculate the kronecker form...\n');
  A,B
  error('sorry')
end;  
%%%%%%

[A,B,Q1a,Z1a]=grsd(A,B);
[A,B,Q1b,Z1b,na]=eigsort(A,B,alpha,option);


Z1=Z1a*Z1b;
Q1=Q1b*Q1a;
nb=n-na;

%
%  converting A           and          B to
%
%      |  A11*inv(B11)   A12 |            |  I   B12 |
%      |   0             A22 |            |  0   B22 |

if na>0,
   Z2=eye(n);Z2(1:na,1:na)=inv(B(1:na,1:na));
   A=A*Z2;
   B=B*Z2;
else
   Z2=eye(n);
end

%  converting A           and          B to
%
%      |  A11*inv(B11)   A12 |            |  I            B12 |
%      |   0              I  |            |  0   inv(A22)*B22 |
% s.t.
%      |eig(A11)*inv(B11)| <= alpha      |eig(inv(A22)*B22)| <1/alpha

if na<n,
   Q2=eye(n);
   Q2(na+1:n,na+1:n)=inv(A(na+1:n,na+1:n));
   A=Q2*A;
   B=Q2*B;
else
   Q2=eye(n);
end

A11=A(1:na,1:na);
B22=B(na+1:n,na+1:n);
A12=A(1:na,na+1:n);
B12=B(1:na,na+1:n);
if na>0 & na<n,   
  if min(abs(eig(B22)))>beta*10
    QQ=lyap(-A11,eye(nb)/B22,A12/B22-A11*B12/B22);
    ZZ = -B12-QQ*B22;
  elseif  min(abs(eig(A11)))>beta*10
    ZZ=lyap(A11\eye(na),-B22,A11\B12-A11\A12*B22);
    QQ = -A12-A11*ZZ;
  else
    error('Both A11 and B22 are singular')
  end
  Q3 = [eye(na) QQ;zeros(nb,na) eye(nb)];
  Z3 = [eye(na) ZZ;zeros(nb,na) eye(nb)];
else   
  Q3=eye(n);Z3=eye(n);
end
A = Q3*A*Z3;
B = Q3*B*Z3;
Q = Q3*Q2*Q1;
Z = Z1*Z2*Z3;

%matrixes opschonen
AA=real(A);BB=real(B);Q=real(Q);Z=real(Z);

for i=1:n,
   for j=1:n,
     if abs(AA(i,j))<beta,AA(i,j)=0;end;
     if abs(BB(i,j))<beta,BB(i,j)=0;end;
   end
end

 
function [AA,BB,Q,Z]=grsd(A,B);

% grsd    Generalized Real Schur Decomposition
%              AA=Q*A*Z
%              BB=Q*B*Z
%         This function transforms a pensil of two real matrixes of size nxn
%         into a upper quasi triangular matrix with 1x1 and 2x2 blocks on 
%         the diagonal, AA, and a upper triagular matrix BB.
%         The ratio between the 1x1 blocks on the diagonal of AA 
%         and the coresponding entry on the diagonal of BB is an eigenvalue 
%         of the pensil A,B
%         The 2x2 blocks on the diagonal of AA correspond with the complex 
%         eigenvalues of the pensil
%
% Syntax:  
%         [AA,BB,Q,Z]=grsd(A,B);
%
% Input:
%  
%  A,B    Two real nxn matrixes.
%  
% Output 
%   AA    QxAxZ upper quasi-triangular real matrix
%               This means that A has 1x1 or 2x2 blocks on the diagonal.
%   BB    QxBxZ upper triangular matrix
%   Q,Z   Transformation matrix
%  


% This function is described in Gene H. Golub, Matrix Computations,
% John Hopkins University Press, second edition,1989
% algorithm 7.7.3 page 403-404;

% L.R.J. Haverkamp 09-01-1995
% copyright (c) 1995 L.R.J. Haverkamp

beta=eps*norm(A)*norm(B);
[A,B,Q,Z]=ghess(A,B);

q=0;
n=length(A);
while q<n,
  % clean subdiagonal
  for i=2:n,
    if abs(A(i,i-1))<beta*(abs(A(i-1,i-1))+abs(A(i,i))),
      A(i,i-1)=0;
    end
  end
  
  % check for upper-quasi-triangularship(determine k)
  k=n; readyflag=0;
  while (readyflag==0 & k>0),
    if k==1
      %'1x1 block in upper left corner'
      k=k-1;
    elseif k>1 & abs(A(k,k-1))<=beta
      %'1x1 block in the middle or down right'
      k=k-1; 
    elseif k>1 &  abs(imag(eig(A(k-1:k,k-1:k),B(k-1:k,k-1:k)))) > beta ,
      %complex eigenvalues?
      if k==2 & abs(A(k,k-1))>= beta 
	%'2x2 block in  upper left corner'
	k=k-2;
      elseif k>2 & abs(A(k,k-1))>=beta & abs(A(k-1,k-2)) < beta 
	%'2x2 block in the middle or down right'
	k=k-2;
      else
	readyflag=1;	    
      end
    else
      readyflag=1;
    end
  end
  
  q=n-k;
  %q is size of largest quasi upper triangular block in lower right corner.
  if q<n,
    % determine p     
    p=k;
    while(p>2 & abs(A(p,p-1))>beta ),p=p-1;end
    if (p==2 & abs(A(p,p-1))>beta ), p=p-1;end
    %    if q==2,
    %    end
    
    A22=A(p:n-q,p:n-q);
    B22=B(p:n-q,p:n-q);
    n22=length(B22);
    zero_flag=0;
    for k=1:n22,  
      if (abs(B22(k,k))<beta),
	zero_flag=1;
	zero_pos=k;
      end
    end
    if zero_flag==1,
      %         fprintf('zerochasing \n')
      [A22,B22,qq,zz]=zerochase(A22,B22,zero_pos);
      QQ=eye(n);QQ(p:n-q,p:n-q)=qq;
      ZZ=eye(n);ZZ(p:n-q,p:n-q)=zz;
      Q=QQ*Q;
      Z=Z*ZZ;
    else
      [A22,B22,qq,zz]=qzstep(A22,B22); 
      QQ=eye(n);QQ(p:n-q,p:n-q)=qq;
      ZZ=eye(n);ZZ(p:n-q,p:n-q)=zz;
    end
    Q=QQ*Q;
    Z=Z*ZZ;
    A=QQ*A*ZZ;
    B=QQ*B*ZZ;
  end
end
AA=A;BB=B;


function [AA,BB,Q,Z,ns] = eigsort(A,B,alpha,option); 
% eigsort In this function, the eigenvalues of the matrix pensil A, 
%         B are sorted, according  the parameter alpha and option. 
%         A and B are expected to be in generalized real schur form, 
%         as obtained from the function grsd. 
%         A and B are transformed by the matrices Q and Z following 
%                      | A11  A12|                 | B11  B12| 
%           Q*A*Z=AA=  | 0    A22| and  Q*B*Z=BB=  | 0    B22| 
%         The eigenvalues on the diagonal of A and B are then 
%         exchanged such that the ones satisfying the condition given 
%         by alpha and option are in the upper left corner, ie A11, 
%         B11, and the rest moved to the lower right corner(A22, B22). 
% 
% Syntax: [AA,BB,Q,Z,ns]=eigsort(A,B,alpha,option) 
% Input: 
%   A,B     Matrix pencil. A and B are square matrices of equal 
%           dimension. A and B are assumed to be in Generalized Real 
%           Schur Form, otherwise eigsort will give bogus output. 
%   alpha   The value of alpha is  used to split the eigenvalues of 
%           the pencil. The default value is alpha=1. 
%   option  Option can have two values; 'circle' or 'halfplane'. 
%           When option equals circle, the eigenvalues 
%           are sorted such that the absolute value of the 
%           eigenvalues in A11 and B22 are smaller than alpha. When 
%           option equals 'halfplane', the eigenvalues are sorted 
%           such that the eigenvalues in A11 and B22 have real parts 
%           less than alpha. 
% 
% Output: 
%   AA,BB   Transformed matrix pensil. 
%   Q,Z     Transformation matrices, satisfying AA=Q*A*Z and BB=Q*B*Z. 
%   ns      Dimension of the matrix As. 
% 
% See also: grsd, ghess, kroneckf. 
 
%  --- This file is eigsort.m generated from the source        --- 
%  ---               eigsort.web, using mweb                   --- 
 
%    Bert Haverkamp 04/05/95 
%    Modified 04/06/98 
%    copyright (c) 1995 Bert Haverkamp 
 
% The algoritme is described in:  P van Dooren, A generalized 
% eigenvalue approach for solving riccati equations 
% SIAM J. Sci. Stat. comput. Vol2 no 2 page 121-135 
% A fortran implementation is described in:  P van Dooren, Algortimm 590 
% ACM transactions on mathematical software, vol 8, no 4, page 376-382 
% the fortran implementation is available at NETLIB as exchqz.f 
 
if nargin>4 
  error('Too many arguments for my little comprehention') 
end 
 
if nargin==4 
  if ~strcmp(option,'halfplane') &  ~strcmp(option,'circle')
    error('unknown value of option given'); 
  end 
end 

if nargin<3 
  option = 'circle';
  alpha = []; 
end 
 
 
if isempty(alpha) 
  if strcmp(option,'circle')
    alpha = 1; 
  else 
    alpha = 0; 
  end 
end 
 
if (size(A,1)~=size(A,2)|size(B,1)~=size(B,2)) 
  error('Matrices A and B have to be square') 
end 
 
if (size(A,1)~=size(B,1)) 
  error('Matrices A and B have to have the same size') 
end 
 
beta = eps * norm(A) * norm(B); 
n = length(A); 
Q = eye(n); 
Z = eye(n); 
sort=[];
nlow=0;
i = 1; 
while i<n, 
  if abs(A(i+1,i))<beta, %[1x1 block] 
    lambda =  eig(A(i,i),B(i,i)); 
    if (strcmp(option,'circle') & (abs(lambda) < alpha)) | ... 
	  (strcmp(option,'halfplane') & (real(lambda) < alpha)) 
      sort = [sort,1]; 
    else 
      sort = [sort,-1]; 
    end 
    i = i+1; 
  else               % [2x2 block]
    lambda =  eig(A(i:i+1,i:i+1),B(i:i+1,i:i+1));
    lambda = lambda(1); 
     if (strcmp(option,'circle') & (abs(lambda) < alpha)) | ... 
	  (strcmp(option,'halfplane') & (real(lambda) < alpha)) 
      sort = [sort,2]; 
    else 
      sort = [sort,-2]; 
    end 
    i = i+2; 
   end 
end 
 if i==n, 
  lambda =  eig(A(i,i),B(i,i)); 
   if (strcmp(option,'circle') & (abs(lambda) < alpha)) | ... 
	 (strcmp(option,'halfplane') & (real(lambda) < alpha)) 
     sort = [sort,1]; 
     nlow = nlow+1; 
   else 
    sort = [sort,-1]; 
   end 
end 
 
  
done_a_swap_flag = 1; 
while done_a_swap_flag==1, 
  done_a_swap_flag = 0; 
  pos = 1; 
  for i = 1:length(sort)-1; 
    if i>1 
      pos = pos+abs(sort(i-1));
    end; 
    if (sort(i)==-1 & sort(i+1)==1), 
      p1 = pos;p2 = pos+1; 
      H = A(p2,p2) * B(p1:p2,p1:p2)-B(p2,p2) * A(p1:p2,p1:p2); 
      z1 = rot90(givens(H(1,1),H(1,2)),3); 
      ZZ = eye(n);ZZ(p1:p2,p1:p2) = z1; 
      B = B * ZZ; 
      q1 = givens(B(p1,p1),B(p2,p1)); 
      QQ = eye(n);QQ(p1:p2,p1:p2) = q1; 
      A = QQ * A * ZZ; 
      B = QQ * B; 
      Q = QQ * Q; 
      Z = Z * ZZ; 
      sort(i) = 1;sort(i+1) = -1; 
      done_a_swap_flag = 1; 
    elseif (sort(i)==-2 & sort(i+1)==1) 
      p1 = pos;p2 = pos+1;p3 = pos+2; 
      H = A(p3,p3) * B(p1:p3,p1:p3)- B(p3,p3) * A(p1:p3,p1:p3); 
      r = givens(H(1,1),H(2,1)); 
      R = eye(3);R(1:2,1:2) = r; 
      H = R * H; 
      z1 = givens(H(2,2),H(2,3)); 
      ZZ1 = eye(n);ZZ1(p2:p3,p2:p3) = rot90(z1,3); 
      H = H * ZZ1(p1:p3,p1:p3); 
      z2 = givens(H(1,1),H(1,2)); 
      ZZ2 = eye(n);ZZ2(p1:p2,p1:p2) = rot90(z2,3); 
      B = B * ZZ1; 
      q1 = givens(B(p2,p2),B(p3,p2)); 
      QQ1 = eye(n);QQ1(p2:p3,p2:p3) = q1; 
      B = QQ1 * B * ZZ2; 
      q2 = givens(B(p1,p1),B(p2,p1)); 
      QQ2 = eye(n);QQ2(p1:p2,p1:p2) = q2; 
      B = QQ2 * B; 
      A = QQ2 * QQ1 * A * ZZ1 * ZZ2; 
      Q = QQ2 * QQ1 * Q; 
      Z = Z * ZZ1 * ZZ2; 
      sort(i) = 1;sort(i+1) = -2; 
      done_a_swap_flag = 1; 
    elseif (sort(i)==-1 & sort(i+1)==2) 
      p1 = n-pos-1;p2 = n-pos;p3 = n-pos+1; 
      A = flipud(fliplr(A')); 
      B = flipud(fliplr(B')); 
      t1 = givens(A(p1,p1),A(p2,p1)); 
      T1 = eye(n);T1(p1:p2,p1:p2) = t1; 
      X = T1 * A; 
      A = T1 * B; 
      B = X; 
      H = A(p3,p3) * B(p1:p3,p1:p3)- B(p3,p3) * A(p1:p3,p1:p3); 
      r = givens(H(1,1),H(2,1)); 
      R = eye(3);R(1:2,1:2) = r; 
      H = R * H; 
      z1 = givens(H(2,2),H(2,3)); 
      ZZ1 = eye(n);ZZ1(p2:p3,p2:p3) = rot90(z1,3); 
      H = H * ZZ1(p1:p3,p1:p3); 
      z2 = givens(H(1,1),H(1,2)); 
      ZZ2 = eye(n);ZZ2(p1:p2,p1:p2) = rot90(z2,3); 
      B = B * ZZ1; 
      q1 = givens(B(p2,p2),B(p3,p2)); 
      QQ1 = eye(n);QQ1(p2:p3,p2:p3) = q1; 
      B = QQ1 * B * ZZ2; 
      q2 = givens(B(p1,p1),B(p2,p1)); 
      QQ2 = eye(n);QQ2(p1:p2,p1:p2) = q2; 
      B = QQ2 * B; 
      A = QQ2 * QQ1 * A * ZZ1 * ZZ2; 
      t2 = givens(A(p3,p2),A(p3,p3)); 
      T2 = eye(n);T2(p2:p3,p2:p3) = rot90(t2,3); 
      X = A * T2; 
      A = B * T2; 
      B = X; 
      A = flipud(fliplr(A')); 
      B = flipud(fliplr(B')); 
      Q = flipud(fliplr((ZZ1 * ZZ2 * T2)')) * Q; 
      Z = Z * flipud(fliplr((QQ2 * QQ1 * T1)')); 
      sort(i) = 2;sort(i+1) = -1; 
      done_a_swap_flag = 1; 
    elseif sort(i)==-2 & sort(i+1)==2, 
      p1 = pos;p2 = pos+1;p3 = pos+2;p4 = pos+3; 
      AMMBMM  =  A(p1,p1)/B(p1,p1); 
      ANMBMM  =  A(p2,p1)/B(p1,p1); 
      AMNBNN  =  A(p1,p2)/B(p2,p2); 
      ANNBNN  =  A(p2,p2)/B(p2,p2); 
      BMNBNN  =  B(p1,p2)/B(p2,p2); 
      A0 = [1;1;1;0]; 
      first_time_flag = 1; 
      while abs(A(p3,p2))>beta | first_time_flag==1, 
	first_time_flag = 0; 
	%QQ1 
	q1 = givens(A0(2),A0(3)); 
	QQ1 = eye(n);QQ1(p2:p3,p2:p3) = q1; 
	A0 = QQ1(p1:p4,p1:p4) * A0; 
	A = QQ1 * A; 
	B = QQ1 * B; 
	%QQ2 
	q2 = givens(A0(1),A0(2)); 
	QQ2 = eye(n);QQ2(p1:p2,p1:p2) = q2; 
	A0 = QQ2(p1:p4,p1:p4) * A0; 
	A = QQ2 * A; 
	B = QQ2 * B; 
	%ZZ1 
	z1 = givens(B(p3,p2),B(p3,p3)); 
	ZZ1 = eye(n);ZZ1(p2:p3,p2:p3) = rot90(z1,3); 
	B = B * ZZ1; 
	A = A * ZZ1; 
	%ZZ2 
	z2 = givens(B(p2,p1),B(p2,p2)); 
	ZZ2 = eye(n);ZZ2(p1:p2,p1:p2) = rot90(z2,3); 
	B = B * ZZ2; 
	A = A * ZZ2; 
	%QQ3 
	q3 = givens(A(p3,p1),A(p4,p1)); 
	QQ3 = eye(n);QQ3(p3:p4,p3:p4) = q3; 
	A = QQ3 * A; 
	B = QQ3 * B; 
	%ZZ3 
	z3 = givens(B(p4,p3),B(p4,p4)); 
	ZZ3 = eye(n);ZZ3(p3:p4,p3:p4) = rot90(z3,3); 
	B = B * ZZ3; 
	A = A * ZZ3; 
	%QQ4 
	q4 = givens(A(p2,p1),A(p3,p1)); 
	QQ4 = eye(n);QQ4(p2:p3,p2:p3) = q4; 
	A = QQ4 * A; 
	B = QQ4 * B; 
	%ZZ4 
	z4 = givens(B(p3,p2),B(p3,p3)); 
	ZZ4 = eye(n);ZZ4(p2:p3,p2:p3) = rot90(z4,3); 
	B = B * ZZ4; 
	A = A * ZZ4; 
	%QQ5 
	q5 = givens(A(p3,p2),A(p4,p2)); 
	QQ5 = eye(n);QQ5(p3:p4,p3:p4) = q5; 
	A = QQ5 * A; 
	B = QQ5 * B; 
	%ZZ5 
	z5 = givens(B(p4,p3),B(p4,p4)); 
	ZZ5 = eye(n);ZZ5(p3:p4,p3:p4) = rot90(z5,3); 
	B = B * ZZ5; 
	A = A * ZZ5; 
	
	A11B11  =  A(p1,p1)/B(p1,p1); 
	A12B22  =  A(p1,p2)/B(p2,p2); 
	A21B11  =  A(p2,p1)/B(p1,p1); 
	A22B22  =  A(p2,p2)/B(p2,p2); 
	B12B22  =  B(p1,p2)/B(p2,p2); 
	x  =  ((AMMBMM-A11B11) * (ANNBNN-A11B11) - AMNBNN * ANMBMM + ... 
	       ANMBMM * BMNBNN * A11B11)/A21B11+A12B22- A11B11 * B12B22; 
	y  =  (A22B22-A11B11)-A21B11 * B12B22-(AMMBMM-A11B11) - ... 
	      (ANNBNN-A11B11) + ANMBMM * BMNBNN; 
	z  =  A(p3,p2)/B(p2,p2); 
	A0 = [x;y;z;0]; 
	Q = QQ5 * QQ4 * QQ3 * QQ2 * QQ1 * Q; 
	Z = Z * ZZ1 * ZZ2 * ZZ3 * ZZ4 * ZZ5; 
      end 
      sort(i) = 2;sort(i+1) = -2; 
      done_a_swap_flag = 1; 
    end 
  end 
end 
ns = sum(sort((sort>0))); 
A = real(A); 
B = real(B); 
Q = real(Q); 
Z = real(Z); 
A = A.*(abs(A)>beta); 
B = B.*(abs(B)>beta); 
AA = A;BB = B; 
 
function [AA,BB,Q,Z]=ghess(A,B);
% function [AA,BB,Q,Z]=ghess(A,B);
% general hessenberg triangular reduction
% input   A  real nxn matrix
%         B  real nxn matrix
% output  AA=Q*A*Z upper hesenberg
%         EE=Q*E*Z upper triangular

% This function is used by grsd.m
% This function uses no other m-files

% This function is described in Gene H. Golub, Matrix Computations,
% John Hopkins University Press, second edition,1989
% algorithm 7.7.1 page 397-399

% L.R.J. Haverkamp 09-01-1995
% copyright (c) 1995 L.R.J. Haverkamp

n=length(A);
Z=eye(n);
[U,R]=qr(B);
B=U'*B;
A=U'*A;
Q=U';
for j=1:n-2,
   for i=n:-1:j+2,
      qi=givens(A(i-1,j),A(i,j));
      Qi=eye(n);Qi(i-1:i,i-1:i)=qi;
      A=Qi*A;
      B=Qi*B;
      Q=Qi*Q;
      zi=givens(B(i,i-1),B(i,i));
      Zi=eye(n);Zi(i-1:i,i-1:i)=rot90(zi,3);
      A=A*Zi;
      B=B*Zi;
      Z=Z*Zi;
   end
end
AA=triu(A);
for i=1:n-1,
   AA(i+1,i)=A(i+1,i);
end
BB=triu(B);


function [AA,BB,Q,Z]=zerochase(A,B,k);
% input:  A  nxn upper hessian
%         B  nxn upper triangular
%         k position of zero on diagonal of B; B(k,k)=0
% output: AA =Q*A*Z, nxn upper hessian, with A(n,n-1)=0
%         BB =Q*B*Z, nxn upper triangular with B(n,n) = zero 
% This function preforms a single zerochase as described in 
% algorithm 7.7.2 of Matrix Computations of G.H.Golub

% this function is used in grsd.m
% this function uses no other m-files


% L.R.J. Haverkamp 11-01-1995
% copyright (c) 1995 L.R.J. Haverkamp

n=length(A);
Q=eye(n);Z=eye(n);
for i=k:n-1, 
   qm=givens(B(i,i+1),B(i+1,i+1));
   Qm=eye(n);Qm(i:i+1,i:i+1)=qm; 
   A=Qm*A;
   B=Qm*B;
   Q=Qm*Q;
   if (i>1), 
      zm=givens(A(i+1,i-1),A(i+1,i));
      Zm=eye(n);Zm(i-1:i,i-1:i)=rot90(zm,3);
      B=B*Zm;
      A=A*Zm;
      Z=Z*Zm;
   end
end
i=n;
zm=givens(A(i,i-1) ,A(i,i));
Zm=eye(n);Zm(i-1:i,i-1:i)=rot90(zm,3);
BB=B*Zm;
AA=A*Zm;
Z=Z*Zm;


function [AA,BB,Q,Z]=qzstep(A,B);
% function [A,B]=qzstep(A,B);
%    AA=Q*A*Z
%    BB=Q*B*Z
% this function preforms a single qz step as described in 
% algorithm 7.7.2 of Matrix Computations of G.H.Golub
% step 2a  a single QZ step
% The routine is slightly adjusted: instead of housholder matrices, givensrotations are used.
% also 1x1 and 2x2 matrices are supported, althought the results are trivial 
%(espesially the 1x1 case :-) )

% this function is used in grsd.m
% this function uses no other m-files

% L.R.J. Haverkamp 11-01-1995
% copyright (c) 1995 L.R.J. Haverkamp

n=length(A);
Q=eye(n);Z=eye(n);
if n==1,
% do nothing ;-)

elseif n==2,
   l=eig(A,B);
   C=A-l(1)*B; 
   Z=givens(C(2,2),C(2,1));
   %Z=[C(2,2) C(2,1);-C(2,1) C(2,2)];
   %Z=Z/norm(Z)
   C*Z;
   AZ= A*Z;
   %Q=[AZ(1,1), AZ(2,1); AZ(2,1), -AZ(1,1)]; 
   %Q=Q/norm(Q);
   Q=givens(AZ(1,1), AZ(2,1));
   %Q=[0 1 ; 1 0];
   %Z=eye(2);
   A=Q*AZ;
   B=Q*B*Z; 
else
   AMMBMM = A(n-1,n-1)/B(n-1,n-1);
   ANMBMM = A(n,n-1)/B(n-1,n-1);
   AMNBNN = A(n-1,n)/B(n,n);
   ANNBNN = A(n,n)/B(n,n);
   BMNBNN = B(n-1,n)/B(n,n);
   A11B11 = A(1,1)/B(1,1);
   A12B22 = A(1,2)/B(2,2);
   A21B11 = A(2,1)/B(1,1);
   A22B22 = A(2,2)/B(2,2);
   B12B22 = B(1,2)/B(2,2);
   x = ((AMMBMM-A11B11)*(ANNBNN-A11B11) - AMNBNN*ANMBMM +...
       ANMBMM*BMNBNN*A11B11)/A21B11+A12B22- A11B11*B12B22;
   y = (A22B22-A11B11)-A21B11*B12B22-(AMMBMM-A11B11) -...
       (ANNBNN-A11B11) + ANMBMM*BMNBNN;
   z = A(3,2)/B(2,2);

   for k=1:n-2,
      A0=[x;y;z];

      %Q1
      q1=givens(A0(2),A0(3));    
      Q1=eye(n);Q1(k+1:k+2,k+1:k+2)=q1;       
      A0=Q1(k:k+2,k:k+2)*A0;  
      A=Q1*A;
      B=Q1*B;

      %Q2
      q2=givens(A0(1),A0(2));
      Q2=eye(n);Q2(k:k+1,k:k+1)=q2; 
      A0=Q2(k:k+2,k:k+2)*A0;
      A=Q2*A;
      B=Q2*B;

      %Z1
      z1=givens(B(k+2,k+1),B(k+2,k+2));
      Z1=eye(n);Z1(k+1:k+2,k+1:k+2)=rot90(z1,3);
      B=B*Z1;
      A=A*Z1;

      %Z2
      z2=givens(B(k+1,k),B(k+1,k+1));
      Z2=eye(n);Z2(k:k+1,k:k+1)=rot90(z2,3);
      B=B*Z2;
      A=A*Z2;
      x=A(k+1,k);
      y=A(k+2,k);
      if (k<n-2), z=A(k+3,k); end;
      Q=Q2*Q1*Q;
      Z=Z*Z1*Z2;
   end
%Q_n-1
      qm=givens(x,y);
      Qm=eye(n);Qm(k+1:k+2,k+1:k+2)=qm; 
      A=Qm*A;
      B=Qm*B;
 %Z_n-1
      zm=givens(B(n,n-1),B(n,n));
      Zm=eye(n);Zm(n-1:n,n-1:n)=rot90(zm,3);
      B=B*Zm;
      A=A*Zm;

Q=Qm*Q;
Z=Z*Zm;
end
AA=A;
BB=B;
