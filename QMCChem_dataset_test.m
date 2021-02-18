Sold = dlmread('Slater_old.dat');
Sold_inv = dlmread('Slater_old_inv.dat');
S = dlmread('Slater.dat');
S_inv = dlmread('Slater_inv.dat');
det(Sold*transpose(Sold_inv))
trace(Sold*transpose(Sold_inv))
norm(Sold*transpose(Sold_inv))

det(S*transpose(S_inv))
trace(S*transpose(S_inv))
norm(S*transpose(S_inv))

