## Algorithm 3 from P. Maponi,
## p. 283, doi:10.1016/j.laa.2006.07.007
clc ## Clear the screen

## Define the matrix to be inverted. This is example 8 from the paper
## In the future this matrix needs to be read from the function call arguments
A=[1,1,-1; ...
   1,1,0; ...
  -1,0,-1];
A0=diag(diag(A));  ## The diagonal part of A

### The modified example that gives all singular updates at some point  
#A=[1,1,1; ...
#   1,1,0; ...
#  -1,0,-1];
#A0=diag(diag(A));  ## The diagonal part of A

### A square uniform distributed random integer matrix with entries in [-1,1]  
#do
#  A=randi([-1,1],3,3);
#  A0=diag(diag(A));  ## The diagonal part of A
#until (det(A)!=0 && det(A0)!=0) ## We need both matrices to be simultaniously non-singular

### A square uniform distributed random float matrix with entries in (0,1)  
#do
#  A=rand(5);
#  A0=diag(diag(A));  ## The diagonal part of A
#until (det(A)!=0 && det(A0)!=0) ## We need both matrices to be simultaniously non-singular

Ar=A-A0;           ## The remainder of A
nCols=columns(Ar); ## The number of coluns of A (M in accompanying PDF)
Id=eye(nCols);
A0inv=eye(nCols);
Ainv=zeros(nCols,nCols);
ylk=zeros(nCols,nCols,nCols);
p=zeros(nCols,1);
breakdown=zeros(nCols,1);

A,A0
printf("Determinant of A  is: %d\n",det(A))
printf("Determinant of A0 is: %d\n",det(A0))

## Calculate the inverse of A0 and populate p-vector
for i=1:nCols
  A0inv(i,i) = 1 / A0(i,i);
  p(i)=i;
endfor

## Calculate all the y0k in M^2 multiplications instead of M^3
for k=1:nCols
  for i=1:nCols
    #printf("(i,k,1) = (%d,%d,1)\n",i,k);
    ylk(i,k,1) = A0inv(i,i) * Ar(i,k);
  endfor
endfor

## Calculate all the ylk from the y0k calculated previously
for l=2:nCols
  ## Calculate break-down conditions and put in a vector
  for j=l-1:nCols
    breakdown(j) = abs(1+ylk(p(j),p(j),l-1));
    #printf("|1 + ylk(%d,%d,%d)|\n", p(j), p(j), l-1);
  endfor
  [val, lbar] = max(breakdown); ## Find the index of the max value
  breakdown=zeros(nCols,1); ## Reset the entries to zero for next l-round
  ## Swap p(l) and p(lbar)
  tmp=p(l-1);
  p(l-1)=p(lbar);
  p(lbar)=tmp;
  for k=l:nCols
    for i=1:nCols
      ylk(i,p(k),l) = ylk(i,p(k),l-1) - (ylk(p(l-1),p(k),l-1)) / (1+ylk(p(l-1),p(l-1),l-1)) * (ylk(i,p(l-1),l-1));
      #printf("ylk(%d,%d,%d) = ylk(%d,%d,%d) - (ylk(%d,%d,%d) / (1+ylk(%d,%d,%d) * (ylk(%d,%d,%d);\n", i,p(k),l,i,p(k),l-1,p(l-1),p(k),l-1,p(l-1),p(l-1),l-1,i,p(l-1),l-1);
    endfor
  endfor
endfor

## Construct A-inverse from A0-inverse and the ylk
Ainv=A0inv;
for l=1:nCols
  Ainv=(Id - ylk(:,p(l),l) * transpose(Id(:,p(l))) / (1 + ylk(p(l),p(l),l))) * Ainv;
  #printf("Ainv=(Id - ylk(:,%d,%d) * transpose(Id(:,%d)) / (1 + ylk(%d,%d,%d))) * Ainv\n",p(l),l,p(l),p(l),p(l),l);
endfor

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
    IdTest
  endif
endif