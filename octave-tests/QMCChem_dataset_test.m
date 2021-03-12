Sold = dlmread('Slater_old.dat');
Sold_inv = dlmread('Slater_old_inv.dat');
S = dlmread('Slater.dat');
S_inv = dlmread('Slater_inv.dat');

printf("\n")
printf("======================================\n")
printf("OLD Slater-matrix and inverse\n")
printf("--------------------------------------\n")
printf("Determinant of S x S_inv : %f\n", det(Sold*Sold_inv))
printf("Trace of S x S_inv       : %f\n", trace(Sold*Sold_inv))
printf("Norm of S x S_inv        : %f\n", norm(Sold*Sold_inv))

printf("\n")
printf("NEW Slater-matrix and inverse\n")
printf("--------------------------------------\n")
printf("Determinant of S x S_inv : %f\n", det(S*S_inv))
printf("Trace of S x S_inv       : %f\n", trace(S*S_inv))
printf("Norm of S x S_inv        : %f\n", norm(S*S_inv))

cutoff = 1e-6;
printf("\n")
printf("Cutoff set to %e: S x S_inv = \n", cutoff)
printf(" -----------------------------------------\n")
dim = columns(S);
res=S*S_inv;
for i = 1:dim
    for j = 1:dim
        if (res(i,j) < cutoff)
            res(i,j) = 0;
        elseif (S(i,j)-1 < cutoff)
            res(i,j) = 1;
        endif
    endfor
endfor
format free;
disp(res);
printf(" =========================================\n")
