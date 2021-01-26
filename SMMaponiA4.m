## Algorithm 4 from P. Maponi,
## p. 283, doi:10.1016/j.laa.2006.07.007
clc ## Clear the screen

% ## Define the matrix to be inverted. This is example 8 from the paper
% ## In the future this matrix needs to be read from the function call arguments
% A=[1,1,-1; ...
%    1,1,0; ...
%   -1,0,-1];

%## The modified example that gives all singular updates at some point
%A=[1,1,1; ...
%   1,1,0; ...
%  -1,0,-1];

## A square uniform distributed random integer matrix with entries in [-1,1]  
do
  A=randi([-5,5],5,5);
until (det(A)!=0) ## We need matrix non-singular

% ## A square uniform distributed random float matrix with entries in (0,1)  
% do
%   A=rand(5);
% until (det(A)!=0) ## We need matrix non-singular


nCols=columns(A); ## The number of coluns of A (M in accompanying PDF)
Id=eye(nCols);
% d=norm(A);
% d=max(eig(A));
d=0.1
A0=d*eye(nCols);
Ar=A-A0;           ## The remainder of A
A0inv=eye(nCols);
Ainv=zeros(nCols,nCols);
ylk=zeros(nCols,nCols,nCols);
breakdown=zeros(nCols,1);
cutOff=1e-10;
p=zeros(nCols,1);
P=zeros(nCols);

## Calculate the inverse of A0 and populate p-vector
for i=1:nCols
  A0inv(i,i) = 1 / A0(i,i);
  p(i)=i;
endfor
A,A0,Ar,A0inv

printf("Determinant of A  is: %d\n",det(A))
printf("Determinant of A0 is: %d\n",det(A0))
printf("Determinant of A0inv is: %d\n",det(A0inv))



## Calculate all the y0k in M^2 multiplications instead of M^3
for k=1:nCols
  for i=1:nCols
    ylk(i,k,1) = A0inv(i,i) * Ar(i,k);
    % printf("ylk(%d,%d,1) = A0inv(%d,%d) * Ar(%d,%d)\n",i,k,i,i,i,k);
  endfor
endfor



## Calculate all the ylk from the y0k calculated previously
for l=2:nCols
  el=Id(:,l-1);
  ## Calculate break-down conditions and put in a vector
  for j=l-1:nCols
    %breakdown(j) = abs( el(j) + ylk(j,l-1,l-1) ) ## Condition from Maponi. Probably not correct!
    l
    ylk(j,j,l-1);
    breakdown(j) = abs( el(j) + ylk(j,j,l-1) )
    % printf("|el(%d) + ylk(%d,%d,%d)|\n", j, j, l-1, l-1);
  endfor
  [val, s] = max(breakdown) ## Find the index of the max value
  % printf("l = %d\ns = %d\n",l-1,s);
  breakdown=zeros(nCols,1); ## Reset the entries to zero for next l-round
  if (s!=l-1) ## Apply partial pivoting
  ## Swap yl-1,k(r) and yl-1,k(s) for all k=l,l+1,...,M
    r=l-1;
    tmp=p(r);
    p(r)=p(s);
    p(s)=tmp;
    for k=l-1:nCols
      tmp=ylk(r,k,l-1);
      ylk(r,k,l-1)=ylk(s,k,l-1);
      % printf("ylk(%d,%d,%d)=ylk(%d,%d,%d)\n",r,k,l-1,s,k,l-1);
      ylk(s,k,l-1)=tmp;
    endfor
    ## Modify yl-1,r and yl-1,s
    er=Id(:,r);
    es=Id(:,s);
    ylk(:,r,l-1) = ylk(:,r,l-1) + es - er;
    ylk(:,s,l-1) = ylk(:,s,l-1) + er - es;
  endif
  ## Compute finally the yl,k
  for k=l:nCols
    for i=1:nCols
      ylk(i,k,l) = ylk(i,k,l-1) - ylk(l-1,k,l-1) / (1 + ylk(l-1,l-1,l-1)) * ylk(i,l-1,l-1);
      % printf("ylk(%d,%d,%d) = ylk(%d,%d,%d) - (ylk(%d,%d,%d) / (1+ylk(%d,%d,%d) * (ylk(%d,%d,%d)\n",i,k,l,i,k,l-1,l-1,k,l-1,l-1,l-1,l-1,i,l-1,l-1);
    endfor
  endfor
endfor



## Construct A-inverse from A0-inverse and the ylk
Ainv=A0inv;
for l=1:nCols
  Ainv=(Id - ylk(:,l,l) * transpose(Id(:,l)) / (1 + ylk(l,l,l))) * Ainv;
  P(:,l)=Id(:,p(l)); ## Construct permutation matrix to swap columns of Ainv
  % printf("Ainv=(Id - ylk(:,%d,%d) * transpose(Id(:,%d)) / (1 + ylk(%d,%d,%d))) * Ainv\n",l,l,l,l,l,l);
endfor
Ainv=transpose(P*transpose(Ainv)); ## Swap the 



## Test if the inverse found is really an inverse (does not work if values are floats)
IdTest=A*Ainv;
if (IdTest==eye(nCols))
   printf("\n");
   printf("Inverse of A^{-1} FOUND!\n");
   Ainv
else
  printf("\n");
  printf("Inverse of A^{-1} NOT found yet.\nRunning another test...\n");
  for i=1:nCols
    for j=1:nCols
      if (abs(IdTest(i,j))<cutOff)
        IdTest(i,j)=0;
      elseif (abs(IdTest(i,j))-1<cutOff)
        IdTest(i,j)=1;
      endif
    endfor
  endfor
  if (IdTest==eye(nCols))
    printf("\n");
    printf("Inverse of A^{-1} FOUND!\n");
    Ainv
  else
    printf("\n");
    printf("Still not found. Giving up!\n");
    Ainv,IdTest
  endif
endif