Sold = dlmread('Slater_old.dat');
Sold_inv = dlmread('Slater_old_inv.dat');
S = dlmread('Slater.dat');
S_inv = dlmread('Slater_inv.dat');

det(Sold*Sold_inv)
trace(Sold*Sold_inv)
norm(Sold*Sold_inv)

det(S*S_inv)
trace(S*S_inv)
norm(S*S_inv)
