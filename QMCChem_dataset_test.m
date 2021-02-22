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
printf("======================================\n")
