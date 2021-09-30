Sold = dlmread('Slater_old.dat');
Sold_inv = dlmread('Slater_old_inv.dat');
S = dlmread('Slater.dat');
S_inv = dlmread('Slater_inv.dat');
dim = columns(Sold);
cutoff = 1e-3;

printf("\n")
printf("============================================\n")
printf("OLD Slater-matrix and inverse\n")
printf("--------------------------------------------\n")
printf("Determinant of S x S_inv - Id : %f\n", det(Sold*Sold_inv-eye(dim)))
printf("Trace of S x S_inv - Id       : %f\n", trace(Sold*Sold_inv-eye(dim)))
printf("Norm of S x S_inv - Id        : %f\n", norm(Sold*Sold_inv-eye(dim)))

printf("\n")
printf("Cutoff set to %e: S x S_inv - Id = \n", cutoff)
printf("--------------------------------------------\n")
res=Sold*Sold_inv-eye(dim);
for i = 1:dim
    for j = 1:dim
        if ( abs(res(i,j)) < cutoff )
            res(i,j) = 0;
        elseif ( abs(res(i,j) - 1) < cutoff )
            res(i,j) = 1;
        endif
    endfor
endfor
format free;
disp(res);
printf("===========================================\n")


printf("\n")
printf("NEW Slater-matrix and inverse\n")
printf("--------------------------------------------\n")
printf("Determinant of S x S_inv - Id : %f\n", det(S*S_inv-eye(dim)))
printf("Trace of S x S_inv - Id       : %f\n", trace(S*S_inv-eye(dim)))
printf("Norm of S x S_inv - Id        : %f\n", norm(S*S_inv-eye(dim)))

printf("\n")
printf("Cutoff set to %e: S x S_inv - Id = \n", cutoff)
printf("--------------------------------------------\n")
res=S*S_inv-eye(dim);
for i = 1:dim
    for j = 1:dim
        if (res(i,j) < cutoff)% && i!=j)
            res(i,j) = 0;
        endif
    endfor
endfor
format free;
disp(res);
printf("===========================================\n")
