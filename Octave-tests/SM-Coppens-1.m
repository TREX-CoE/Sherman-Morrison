## Based on Algorithm 3 from P. Maponi,
## p. 283, doi:10.1016/j.laa.2006.07.007

clc ## Clear the screen
## Define the matrix to be inverted. This is example 8 from the paper
## In the future this matrix needs to be read from the function call arguments
#A=[1,1,-1; ...
#   1,1,0; ...
#  -1,0,-1];
#A0=diag(diag(A));  ## The diagonal part of A

### The modified example that gives all singular updates at some point  
#A=[1,1,1; ...
#   1,1,0; ...
#  -1,0,-1];
#A0=diag(diag(A));  ## The diagonal part of A

## A uniform distributed random square (5x5) integer matrix with entries in [-1,1]  
do
  A=randi([-1,1],5,5);
  #A=rand(5);
  A0=diag(diag(A));  ## The diagonal part of A
until (det(A)!=0 && det(A0)!=0) ## We need both matrices to be simultaniously non-singular
A
printf("Determinant of A  is: %d\n",det(A))
printf("Determinant of A0 is: %d\n\n",det(A0))

Ar=A-A0;           ## The remainder of A
nCols=columns(Ar); ## The number of coluns of A (M in accompanying PDF)
Id=eye(nCols);     ## Identity matrix, used for the Cartesian basis vectors
Ainv=zeros(nCols,nCols,nCols+1);  ## 3d matrix to store intermediate inverse of A
U=zeros(nCols,nCols,nCols);       ## 3d matrix to store the nCols rank-1 updates
a=zeros(1,nCols);     ## Vector containing the break-down values
used=zeros(1,nCols);  ## Vector to keep track of which updates have been used

## Loop to calculate A0_inv and the rank-1 updates and store in U
Ainv(:,:,1)=eye(nCols);
for i=1:nCols
  Ainv(i,i,1)=1/A0(i,i);
  U(:,:,i)=Ar(:,i)*transpose(Id(:,i));
endfor

## Here starts the calculation of the inverse of A
for i=1:nCols ## Outer loop iterated over the intermediates A_l^-1, M times in total

  ## Here the break-down values are calculated for each intermediate inverse
  for j=1:nCols
    a(j)=1 + transpose(Id(:,j)) * Ainv(:,:,i) * Ar(:,j); ## First time all elmts 1 because A0_inv is diagonal
    printf("a(%d) = %f\n",j,a(j))
    ## Select next update ONLY if break-down value != 0 AND not yet used
    if (a(j)!=0 && used(j)!=1);
      k=j;
      break;
    endif
  endfor

  ## Do the actual S-M update
  Ainv(:,:,i+1) = Ainv(:,:,i) - Ainv(:,:,i) * U(:,:,k) * Ainv(:,:,i) / a(k);
  used(k)=1;  ## Mark this update as used 
endfor

## Test if the inverse found is really an inverse (does not work if values are floats)
if (A*Ainv(:,:,nCols+1)==eye(nCols))
  printf("\n");
  printf("Inverse found. A^{-1}:\n");
  Ainv(:,:,nCols+1)
else
  printf("\n");
  printf("Inverse NOT found! A*A^{-1}:\n");
  A*Ainv(:,:,nCols+1)
endif
