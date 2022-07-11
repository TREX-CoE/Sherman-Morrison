#! /bin/octave -qf

data_anthony=load('ANTHONY.dat');
data_naive=load('NAIVE.dat');
data_later=load('LATER.dat');
data_split=load('SPLITTING.dat');
data_blocked=load('BLOCKED.dat');
data_lapack=load('MKL_LAPACK.dat');
data_wb2=load('WB2.dat');
data_wb3=load('WB3.dat');

indcs=(data_anthony(:,5)==0); % select cycles that passed
anthony_pass_all=data_anthony(indcs,:);
dlmwrite('anthony_pass_all.dat',anthony_pass_all, ' ')

indcs=(data_anthony(:,5)!=0); % select cycles that failed
anthony_fail_all=data_anthony(indcs,:);
dlmwrite('anthony_fail_all.dat',anthony_fail_all, ' ')

indcs=(anthony_pass_all(:,2)==1); % select cycles that passed containing 1 upd
anthony_pass_1=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_1.dat',anthony_pass_1, ' ')
indcs=(anthony_fail_all(:,2)==1); % select cycles that failed containing 1 upd
anthony_fail_1=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_1.dat',anthony_fail_1, ' ')

indcs=(anthony_pass_all(:,2)==2); % select cycles that passed containing 2 upd
anthony_pass_2=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_2.dat',anthony_pass_2, ' ')
indcs=(anthony_fail_all(:,2)==2); % select cycles that failed containing 2 upd
anthony_fail_2=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_2.dat',anthony_fail_2, ' ')

indcs=(anthony_pass_all(:,2)==3); % select cycles that passed containing 3 upd
anthony_pass_3=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_3.dat',anthony_pass_3, ' ')
indcs=(anthony_fail_all(:,2)==3); % select cycles that failed containing 3 upd
anthony_fail_3=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_3.dat',anthony_fail_3, ' ')

indcs=(anthony_pass_all(:,2)==4); % select cycles that passed containing 4 upd
anthony_pass_4=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_4.dat',anthony_pass_4, ' ')
indcs=(anthony_fail_all(:,2)==4); % select cycles that failed containing 4 upd
anthony_fail_4=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_4.dat',anthony_fail_4, ' ')

indcs=(anthony_pass_all(:,2)==5); % select cycles that passed containing 5 upd
anthony_pass_5=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_5.dat',anthony_pass_5, ' ')
indcs=(anthony_fail_all(:,2)==5); % select cycles that failed containing 5 upd
anthony_fail_5=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_5.dat',anthony_fail_5, ' ')

indcs=(anthony_pass_all(:,2)==6); % select cycles that passed containing 6 upd
anthony_pass_6=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_6.dat',anthony_pass_6, ' ')
indcs=(anthony_fail_all(:,2)==6); % select cycles that failed containing 6 upd
anthony_fail_6=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_6.dat',anthony_fail_6, ' ')

indcs=(anthony_pass_all(:,2)==7); % select cycles that passed containing 7 upd
anthony_pass_7=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_7.dat',anthony_pass_7, ' ')
indcs=(anthony_fail_all(:,2)==7); % select cycles that failed containing 7 upd
anthony_fail_7=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_7.dat',anthony_fail_7, ' ')

indcs=(anthony_pass_all(:,2)==8); % select cycles that passed containing 8 upd
anthony_pass_8=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_8.dat',anthony_pass_8, ' ')
indcs=(anthony_fail_all(:,2)==8); % select cycles that failed containing 8 upd
anthony_fail_8=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_8.dat',anthony_fail_8, ' ')

indcs=(anthony_pass_all(:,2)==9); % select cycles that passed containing 9 upd
anthony_pass_9=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_9.dat',anthony_pass_9, ' ')
indcs=(anthony_fail_all(:,2)==9); % select cycles that failed containing 9 upd
anthony_fail_9=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_9.dat',anthony_fail_9, ' ')

indcs=(anthony_pass_all(:,2)==10); % select cycles that passed containing 10 upd
anthony_pass_10=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_10.dat',anthony_pass_10, ' ')
indcs=(anthony_fail_all(:,2)==10); % select cycles that failed containing 10 upd
anthony_fail_10=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_10.dat',anthony_fail_10, ' ')

indcs=(anthony_pass_all(:,2)==11); % select cycles that passed containing 11 upd
anthony_pass_11=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_11.dat',anthony_pass_11, ' ')
indcs=(anthony_fail_all(:,2)==11); % select cycles that failed containing 11 upd
anthony_fail_11=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_11.dat',anthony_fail_11, ' ')

indcs=(anthony_pass_all(:,2)==12); % select cycles that passed containing 12 upd
anthony_pass_12=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_12.dat',anthony_pass_12, ' ')
indcs=(anthony_fail_all(:,2)==12); % select cycles that failed containing 12 upd
anthony_fail_12=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_12.dat',anthony_fail_12, ' ')

indcs=(anthony_pass_all(:,2)==13); % select cycles that passed containing 13 upd
anthony_pass_13=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_13.dat',anthony_pass_13, ' ')
indcs=(anthony_fail_all(:,2)==13); % select cycles that failed containing 13 upd
anthony_fail_13=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_13.dat',anthony_fail_13, ' ')

indcs=(anthony_pass_all(:,2)==14); % select cycles that passed containing 14 upd
anthony_pass_14=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_14.dat',anthony_pass_14, ' ')
indcs=(anthony_fail_all(:,2)==14); % select cycles that failed containing 14 upd
anthony_fail_14=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_14.dat',anthony_fail_14, ' ')

indcs=(anthony_pass_all(:,2)==15); % select cycles that passed containing 15 upd
anthony_pass_15=anthony_pass_all(indcs,:);
dlmwrite('anthony_pass_15.dat',anthony_pass_15, ' ')
indcs=(anthony_fail_all(:,2)==15); % select cycles that failed containing 15 upd
anthony_fail_15=anthony_fail_all(indcs,:);
dlmwrite('anthony_fail_15.dat',anthony_fail_15, ' ')


indcs=(data_naive(:,5)==0); % select cycles that passed
naive_pass_all=data_naive(indcs,:);
dlmwrite('naive_pass_all.dat',naive_pass_all, ' ')

indcs=(data_naive(:,5)!=0); % select cycles that failed
naive_fail_all=data_naive(indcs,:);
dlmwrite('naive_fail_all.dat',naive_fail_all, ' ')

indcs=(naive_pass_all(:,2)==1); % select cycles that passed containing 1 upd
naive_pass_1=naive_pass_all(indcs,:);
dlmwrite('naive_pass_1.dat',naive_pass_1, ' ')
indcs=(naive_fail_all(:,2)==1); % select cycles that failed containing 1 upd
naive_fail_1=naive_fail_all(indcs,:);
dlmwrite('naive_fail_1.dat',naive_fail_1, ' ')

indcs=(naive_pass_all(:,2)==2); % select cycles that passed containing 2 upd
naive_pass_2=naive_pass_all(indcs,:);
dlmwrite('naive_pass_2.dat',naive_pass_2, ' ')
indcs=(naive_fail_all(:,2)==2); % select cycles that failed containing 2 upd
naive_fail_2=naive_fail_all(indcs,:);
dlmwrite('naive_fail_2.dat',naive_fail_2, ' ')

indcs=(naive_pass_all(:,2)==3); % select cycles that passed containing 3 upd
naive_pass_3=naive_pass_all(indcs,:);
dlmwrite('naive_pass_3.dat',naive_pass_3, ' ')
indcs=(naive_fail_all(:,2)==3); % select cycles that failed containing 3 upd
naive_fail_3=naive_fail_all(indcs,:);
dlmwrite('naive_fail_3.dat',naive_fail_3, ' ')

indcs=(naive_pass_all(:,2)==4); % select cycles that passed containing 4 upd
naive_pass_4=naive_pass_all(indcs,:);
dlmwrite('naive_pass_4.dat',naive_pass_4, ' ')
indcs=(naive_fail_all(:,2)==4); % select cycles that failed containing 4 upd
naive_fail_4=naive_fail_all(indcs,:);
dlmwrite('naive_fail_4.dat',naive_fail_4, ' ')

indcs=(naive_pass_all(:,2)==5); % select cycles that passed containing 5 upd
naive_pass_5=naive_pass_all(indcs,:);
dlmwrite('naive_pass_5.dat',naive_pass_5, ' ')
indcs=(naive_fail_all(:,2)==5); % select cycles that failed containing 5 upd
naive_fail_5=naive_fail_all(indcs,:);
dlmwrite('naive_fail_5.dat',naive_fail_5, ' ')

indcs=(naive_pass_all(:,2)==6); % select cycles that passed containing 6 upd
naive_pass_6=naive_pass_all(indcs,:);
dlmwrite('naive_pass_6.dat',naive_pass_6, ' ')
indcs=(naive_fail_all(:,2)==6); % select cycles that failed containing 6 upd
naive_fail_6=naive_fail_all(indcs,:);
dlmwrite('naive_fail_6.dat',naive_fail_6, ' ')

indcs=(naive_pass_all(:,2)==7); % select cycles that passed containing 7 upd
naive_pass_7=naive_pass_all(indcs,:);
dlmwrite('naive_pass_7.dat',naive_pass_7, ' ')
indcs=(naive_fail_all(:,2)==7); % select cycles that failed containing 7 upd
naive_fail_7=naive_fail_all(indcs,:);
dlmwrite('naive_fail_7.dat',naive_fail_7, ' ')

indcs=(naive_pass_all(:,2)==8); % select cycles that passed containing 8 upd
naive_pass_8=naive_pass_all(indcs,:);
dlmwrite('naive_pass_8.dat',naive_pass_8, ' ')
indcs=(naive_fail_all(:,2)==8); % select cycles that failed containing 8 upd
naive_fail_8=naive_fail_all(indcs,:);
dlmwrite('naive_fail_8.dat',naive_fail_8, ' ')

indcs=(naive_pass_all(:,2)==9); % select cycles that passed containing 9 upd
naive_pass_9=naive_pass_all(indcs,:);
dlmwrite('naive_pass_9.dat',naive_pass_9, ' ')
indcs=(naive_fail_all(:,2)==9); % select cycles that failed containing 9 upd
naive_fail_9=naive_fail_all(indcs,:);
dlmwrite('naive_fail_9.dat',naive_fail_9, ' ')

indcs=(naive_pass_all(:,2)==10); % select cycles that passed containing 10 upd
naive_pass_10=naive_pass_all(indcs,:);
dlmwrite('naive_pass_10.dat',naive_pass_10, ' ')
indcs=(naive_fail_all(:,2)==10); % select cycles that failed containing 10 upd
naive_fail_10=naive_fail_all(indcs,:);
dlmwrite('naive_fail_10.dat',naive_fail_10, ' ')

indcs=(naive_pass_all(:,2)==11); % select cycles that passed containing 11 upd
naive_pass_11=naive_pass_all(indcs,:);
dlmwrite('naive_pass_11.dat',naive_pass_11, ' ')
indcs=(naive_fail_all(:,2)==11); % select cycles that failed containing 11 upd
naive_fail_11=naive_fail_all(indcs,:);
dlmwrite('naive_fail_11.dat',naive_fail_11, ' ')

indcs=(naive_pass_all(:,2)==12); % select cycles that passed containing 12 upd
naive_pass_12=naive_pass_all(indcs,:);
dlmwrite('naive_pass_12.dat',naive_pass_12, ' ')
indcs=(naive_fail_all(:,2)==12); % select cycles that failed containing 12 upd
naive_fail_12=naive_fail_all(indcs,:);
dlmwrite('naive_fail_12.dat',naive_fail_12, ' ')

indcs=(naive_pass_all(:,2)==13); % select cycles that passed containing 13 upd
naive_pass_13=naive_pass_all(indcs,:);
dlmwrite('naive_pass_13.dat',naive_pass_13, ' ')
indcs=(naive_fail_all(:,2)==13); % select cycles that failed containing 13 upd
naive_fail_13=naive_fail_all(indcs,:);
dlmwrite('naive_fail_13.dat',naive_fail_13, ' ')

indcs=(naive_pass_all(:,2)==14); % select cycles that passed containing 14 upd
naive_pass_14=naive_pass_all(indcs,:);
dlmwrite('naive_pass_14.dat',naive_pass_14, ' ')
indcs=(naive_fail_all(:,2)==14); % select cycles that failed containing 14 upd
naive_fail_14=naive_fail_all(indcs,:);
dlmwrite('naive_fail_14.dat',naive_fail_14, ' ')

indcs=(naive_pass_all(:,2)==15); % select cycles that passed containing 15 upd
naive_pass_15=naive_pass_all(indcs,:);
dlmwrite('naive_pass_15.dat',naive_pass_15, ' ')
indcs=(naive_fail_all(:,2)==15); % select cycles that failed containing 15 upd
naive_fail_15=naive_fail_all(indcs,:);
dlmwrite('naive_fail_15.dat',naive_fail_15, ' ')

indcs=(data_later(:,5)==0); % select cycles that passed
later_pass_all=data_later(indcs,:);
dlmwrite('later_pass_all.dat',later_pass_all, ' ')

indcs=(data_later(:,5)!=0); % select cycles that failed
later_fail_all=data_later(indcs,:);
dlmwrite('later_fail_all.dat',later_fail_all, ' ')

indcs=(later_pass_all(:,2)==1); % select cycles that passed containing 1 upd
later_pass_1=later_pass_all(indcs,:);
dlmwrite('later_pass_1.dat',later_pass_1, ' ')
indcs=(later_fail_all(:,2)==1); % select cycles that failed containing 1 upd
later_fail_1=later_fail_all(indcs,:);
dlmwrite('later_fail_1.dat',later_fail_1, ' ')

indcs=(later_pass_all(:,2)==2); % select cycles that passed containing 2 upd
later_pass_2=later_pass_all(indcs,:);
dlmwrite('later_pass_2.dat',later_pass_2, ' ')
indcs=(later_fail_all(:,2)==2); % select cycles that failed containing 2 upd
later_fail_2=later_fail_all(indcs,:);
dlmwrite('later_fail_2.dat',later_fail_2, ' ')

indcs=(later_pass_all(:,2)==3); % select cycles that passed containing 3 upd
later_pass_3=later_pass_all(indcs,:);
dlmwrite('later_pass_3.dat',later_pass_3, ' ')
indcs=(later_fail_all(:,2)==3); % select cycles that failed containing 3 upd
later_fail_3=later_fail_all(indcs,:);
dlmwrite('later_fail_3.dat',later_fail_3, ' ')

indcs=(later_pass_all(:,2)==4); % select cycles that passed containing 4 upd
later_pass_4=later_pass_all(indcs,:);
dlmwrite('later_pass_4.dat',later_pass_4, ' ')
indcs=(later_fail_all(:,2)==4); % select cycles that failed containing 4 upd
later_fail_4=later_fail_all(indcs,:);
dlmwrite('later_fail_4.dat',later_fail_4, ' ')

indcs=(later_pass_all(:,2)==5); % select cycles that passed containing 5 upd
later_pass_5=later_pass_all(indcs,:);
dlmwrite('later_pass_5.dat',later_pass_5, ' ')
indcs=(later_fail_all(:,2)==5); % select cycles that failed containing 5 upd
later_fail_5=later_fail_all(indcs,:);
dlmwrite('later_fail_5.dat',later_fail_5, ' ')

indcs=(later_pass_all(:,2)==6); % select cycles that passed containing 6 upd
later_pass_6=later_pass_all(indcs,:);
dlmwrite('later_pass_6.dat',later_pass_6, ' ')
indcs=(later_fail_all(:,2)==6); % select cycles that failed containing 6 upd
later_fail_6=later_fail_all(indcs,:);
dlmwrite('later_fail_6.dat',later_fail_6, ' ')

indcs=(later_pass_all(:,2)==7); % select cycles that passed containing 7 upd
later_pass_7=later_pass_all(indcs,:);
dlmwrite('later_pass_7.dat',later_pass_7, ' ')
indcs=(later_fail_all(:,2)==7); % select cycles that failed containing 7 upd
later_fail_7=later_fail_all(indcs,:);
dlmwrite('later_fail_7.dat',later_fail_7, ' ')

indcs=(later_pass_all(:,2)==8); % select cycles that passed containing 8 upd
later_pass_8=later_pass_all(indcs,:);
dlmwrite('later_pass_8.dat',later_pass_8, ' ')
indcs=(later_fail_all(:,2)==8); % select cycles that failed containing 8 upd
later_fail_8=later_fail_all(indcs,:);
dlmwrite('later_fail_8.dat',later_fail_8, ' ')

indcs=(later_pass_all(:,2)==9); % select cycles that passed containing 9 upd
later_pass_9=later_pass_all(indcs,:);
dlmwrite('later_pass_9.dat',later_pass_9, ' ')
indcs=(later_fail_all(:,2)==9); % select cycles that failed containing 9 upd
later_fail_9=later_fail_all(indcs,:);
dlmwrite('later_fail_9.dat',later_fail_9, ' ')

indcs=(later_pass_all(:,2)==10); % select cycles that passed containing 10 upd
later_pass_10=later_pass_all(indcs,:);
dlmwrite('later_pass_10.dat',later_pass_10, ' ')
indcs=(later_fail_all(:,2)==10); % select cycles that failed containing 10 upd
later_fail_10=later_fail_all(indcs,:);
dlmwrite('later_fail_10.dat',later_fail_10, ' ')

indcs=(later_pass_all(:,2)==11); % select cycles that passed containing 11 upd
later_pass_11=later_pass_all(indcs,:);
dlmwrite('later_pass_11.dat',later_pass_11, ' ')
indcs=(later_fail_all(:,2)==11); % select cycles that failed containing 11 upd
later_fail_11=later_fail_all(indcs,:);
dlmwrite('later_fail_11.dat',later_fail_11, ' ')

indcs=(later_pass_all(:,2)==12); % select cycles that passed containing 12 upd
later_pass_12=later_pass_all(indcs,:);
dlmwrite('later_pass_12.dat',later_pass_12, ' ')
indcs=(later_fail_all(:,2)==12); % select cycles that failed containing 12 upd
later_fail_12=later_fail_all(indcs,:);
dlmwrite('later_fail_12.dat',later_fail_12, ' ')

indcs=(later_pass_all(:,2)==13); % select cycles that passed containing 13 upd
later_pass_13=later_pass_all(indcs,:);
dlmwrite('later_pass_13.dat',later_pass_13, ' ')
indcs=(later_fail_all(:,2)==13); % select cycles that failed containing 13 upd
later_fail_13=later_fail_all(indcs,:);
dlmwrite('later_fail_13.dat',later_fail_13, ' ')

indcs=(later_pass_all(:,2)==14); % select cycles that passed containing 14 upd
later_pass_14=later_pass_all(indcs,:);
dlmwrite('later_pass_14.dat',later_pass_14, ' ')
indcs=(later_fail_all(:,2)==14); % select cycles that failed containing 14 upd
later_fail_14=later_fail_all(indcs,:);
dlmwrite('later_fail_14.dat',later_fail_14, ' ')

indcs=(later_pass_all(:,2)==15); % select cycles that passed containing 15 upd
later_pass_15=later_pass_all(indcs,:);
dlmwrite('later_pass_15.dat',later_pass_15, ' ')
indcs=(later_fail_all(:,2)==15); % select cycles that failed containing 15 upd
later_fail_15=later_fail_all(indcs,:);
dlmwrite('later_fail_15.dat',later_fail_15, ' ')


indcs=(data_split(:,5)==0); % select cycles that passed
split_pass_all=data_split(indcs,:);
dlmwrite('split_pass_all.dat',split_pass_all, ' ')

indcs=(data_split(:,5)!=0); % select cycles that failed
split_fail_all=data_split(indcs,:);
dlmwrite('split_fail_all.dat',split_fail_all, ' ')

indcs=(split_pass_all(:,2)==1); % select cycles that passed containing 1 upd
split_pass_1=split_pass_all(indcs,:);
dlmwrite('split_pass_1.dat',split_pass_1, ' ')
indcs=(split_fail_all(:,2)==1); % select cycles that failed containing 1 upd
split_fail_1=split_fail_all(indcs,:);
dlmwrite('split_fail_1.dat',split_fail_1, ' ')

indcs=(split_pass_all(:,2)==2); % select cycles that passed containing 2 upd
split_pass_2=split_pass_all(indcs,:);
dlmwrite('split_pass_2.dat',split_pass_2, ' ')
indcs=(split_fail_all(:,2)==2); % select cycles that failed containing 2 upd
split_fail_2=split_fail_all(indcs,:);
dlmwrite('split_fail_2.dat',split_fail_2, ' ')

indcs=(split_pass_all(:,2)==3); % select cycles that passed containing 3 upd
split_pass_3=split_pass_all(indcs,:);
dlmwrite('split_pass_3.dat',split_pass_3, ' ')
indcs=(split_fail_all(:,2)==3); % select cycles that failed containing 3 upd
split_fail_3=split_fail_all(indcs,:);
dlmwrite('split_fail_3.dat',split_fail_3, ' ')

indcs=(split_pass_all(:,2)==4); % select cycles that passed containing 4 upd
split_pass_4=split_pass_all(indcs,:);
dlmwrite('split_pass_4.dat',split_pass_4, ' ')
indcs=(split_fail_all(:,2)==4); % select cycles that failed containing 4 upd
split_fail_4=split_fail_all(indcs,:);
dlmwrite('split_fail_4.dat',split_fail_4, ' ')

indcs=(split_pass_all(:,2)==5); % select cycles that passed containing 5 upd
split_pass_5=split_pass_all(indcs,:);
dlmwrite('split_pass_5.dat',split_pass_5, ' ')
indcs=(split_fail_all(:,2)==5); % select cycles that failed containing 5 upd
split_fail_5=split_fail_all(indcs,:);
dlmwrite('split_fail_5.dat',split_fail_5, ' ')

indcs=(split_pass_all(:,2)==6); % select cycles that passed containing 6 upd
split_pass_6=split_pass_all(indcs,:);
dlmwrite('split_pass_6.dat',split_pass_6, ' ')
indcs=(split_fail_all(:,2)==6); % select cycles that failed containing 6 upd
split_fail_6=split_fail_all(indcs,:);
dlmwrite('split_fail_6.dat',split_fail_6, ' ')

indcs=(split_pass_all(:,2)==7); % select cycles that passed containing 7 upd
split_pass_7=split_pass_all(indcs,:);
dlmwrite('split_pass_7.dat',split_pass_7, ' ')
indcs=(split_fail_all(:,2)==7); % select cycles that failed containing 7 upd
split_fail_7=split_fail_all(indcs,:);
dlmwrite('split_fail_7.dat',split_fail_7, ' ')

indcs=(split_pass_all(:,2)==8); % select cycles that passed containing 8 upd
split_pass_8=split_pass_all(indcs,:);
dlmwrite('split_pass_8.dat',split_pass_8, ' ')
indcs=(split_fail_all(:,2)==8); % select cycles that failed containing 8 upd
split_fail_8=split_fail_all(indcs,:);
dlmwrite('split_fail_8.dat',split_fail_8, ' ')

indcs=(split_pass_all(:,2)==9); % select cycles that passed containing 9 upd
split_pass_9=split_pass_all(indcs,:);
dlmwrite('split_pass_9.dat',split_pass_9, ' ')
indcs=(split_fail_all(:,2)==9); % select cycles that failed containing 9 upd
split_fail_9=split_fail_all(indcs,:);
dlmwrite('split_fail_9.dat',split_fail_9, ' ')

indcs=(split_pass_all(:,2)==10); % select cycles that passed containing 10 upd
split_pass_10=split_pass_all(indcs,:);
dlmwrite('split_pass_10.dat',split_pass_10, ' ')
indcs=(split_fail_all(:,2)==10); % select cycles that failed containing 10 upd
split_fail_10=split_fail_all(indcs,:);
dlmwrite('split_fail_10.dat',split_fail_10, ' ')

indcs=(split_pass_all(:,2)==11); % select cycles that passed containing 11 upd
split_pass_11=split_pass_all(indcs,:);
dlmwrite('split_pass_11.dat',split_pass_11, ' ')
indcs=(split_fail_all(:,2)==11); % select cycles that failed containing 11 upd
split_fail_11=split_fail_all(indcs,:);
dlmwrite('split_fail_11.dat',split_fail_11, ' ')

indcs=(split_pass_all(:,2)==12); % select cycles that passed containing 12 upd
split_pass_12=split_pass_all(indcs,:);
dlmwrite('split_pass_12.dat',split_pass_12, ' ')
indcs=(split_fail_all(:,2)==12); % select cycles that failed containing 12 upd
split_fail_12=split_fail_all(indcs,:);
dlmwrite('split_fail_12.dat',split_fail_12, ' ')

indcs=(split_pass_all(:,2)==13); % select cycles that passed containing 13 upd
split_pass_13=split_pass_all(indcs,:);
dlmwrite('split_pass_13.dat',split_pass_13, ' ')
indcs=(split_fail_all(:,2)==13); % select cycles that failed containing 13 upd
split_fail_13=split_fail_all(indcs,:);
dlmwrite('split_fail_13.dat',split_fail_13, ' ')

indcs=(split_pass_all(:,2)==14); % select cycles that passed containing 14 upd
split_pass_14=split_pass_all(indcs,:);
dlmwrite('split_pass_14.dat',split_pass_14, ' ')
indcs=(split_fail_all(:,2)==14); % select cycles that failed containing 14 upd
split_fail_14=split_fail_all(indcs,:);
dlmwrite('split_fail_14.dat',split_fail_14, ' ')

indcs=(split_pass_all(:,2)==15); % select cycles that passed containing 15 upd
split_pass_15=split_pass_all(indcs,:);
dlmwrite('split_pass_15.dat',split_pass_15, ' ')
indcs=(split_fail_all(:,2)==15); % select cycles that failed containing 15 upd
split_fail_15=split_fail_all(indcs,:);
dlmwrite('split_fail_15.dat',split_fail_15, ' ')


indcs=(data_blocked(:,5)==0); % select cycles that passed
blocked_pass_all=data_blocked(indcs,:);
dlmwrite('blocked_pass_all.dat',blocked_pass_all, ' ')

indcs=(data_blocked(:,5)!=0); % select cycles that failed
blocked_fail_all=data_blocked(indcs,:);
dlmwrite('blocked_fail_all.dat',blocked_fail_all, ' ')

indcs=(blocked_pass_all(:,2)==1); % select cycles that passed containing 1 upd
blocked_pass_1=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_1.dat',blocked_pass_1, ' ')
indcs=(blocked_fail_all(:,2)==1); % select cycles that failed containing 1 upd
blocked_fail_1=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_1.dat',blocked_fail_1, ' ')

indcs=(blocked_pass_all(:,2)==2); % select cycles that passed containing 2 upd
blocked_pass_2=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_2.dat',blocked_pass_2, ' ')
indcs=(blocked_fail_all(:,2)==2); % select cycles that failed containing 2 upd
blocked_fail_2=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_2.dat',blocked_fail_2, ' ')

indcs=(blocked_pass_all(:,2)==3); % select cycles that passed containing 3 upd
blocked_pass_3=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_3.dat',blocked_pass_3, ' ')
indcs=(blocked_fail_all(:,2)==3); % select cycles that failed containing 3 upd
blocked_fail_3=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_3.dat',blocked_fail_3, ' ')

indcs=(blocked_pass_all(:,2)==4); % select cycles that passed containing 4 upd
blocked_pass_4=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_4.dat',blocked_pass_4, ' ')
indcs=(blocked_fail_all(:,2)==4); % select cycles that failed containing 4 upd
blocked_fail_4=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_4.dat',blocked_fail_4, ' ')

indcs=(blocked_pass_all(:,2)==5); % select cycles that passed containing 5 upd
blocked_pass_5=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_5.dat',blocked_pass_5, ' ')
indcs=(blocked_fail_all(:,2)==5); % select cycles that failed containing 5 upd
blocked_fail_5=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_5.dat',blocked_fail_5, ' ')

indcs=(blocked_pass_all(:,2)==6); % select cycles that passed containing 6 upd
blocked_pass_6=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_6.dat',blocked_pass_6, ' ')
indcs=(blocked_fail_all(:,2)==6); % select cycles that failed containing 6 upd
blocked_fail_6=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_6.dat',blocked_fail_6, ' ')

indcs=(blocked_pass_all(:,2)==7); % select cycles that passed containing 7 upd
blocked_pass_7=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_7.dat',blocked_pass_7, ' ')
indcs=(blocked_fail_all(:,2)==7); % select cycles that failed containing 7 upd
blocked_fail_7=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_7.dat',blocked_fail_7, ' ')

indcs=(blocked_pass_all(:,2)==8); % select cycles that passed containing 8 upd
blocked_pass_8=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_8.dat',blocked_pass_8, ' ')
indcs=(blocked_fail_all(:,2)==8); % select cycles that failed containing 8 upd
blocked_fail_8=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_8.dat',blocked_fail_8, ' ')

indcs=(blocked_pass_all(:,2)==9); % select cycles that passed containing 9 upd
blocked_pass_9=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_9.dat',blocked_pass_9, ' ')
indcs=(blocked_fail_all(:,2)==9); % select cycles that failed containing 9 upd
blocked_fail_9=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_9.dat',blocked_fail_9, ' ')

indcs=(blocked_pass_all(:,2)==10); % select cycles that passed containing 10 upd
blocked_pass_10=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_10.dat',blocked_pass_10, ' ')
indcs=(blocked_fail_all(:,2)==10); % select cycles that failed containing 10 upd
blocked_fail_10=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_10.dat',blocked_fail_10, ' ')

indcs=(blocked_pass_all(:,2)==11); % select cycles that passed containing 11 upd
blocked_pass_11=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_11.dat',blocked_pass_11, ' ')
indcs=(blocked_fail_all(:,2)==11); % select cycles that failed containing 11 upd
blocked_fail_11=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_11.dat',blocked_fail_11, ' ')

indcs=(blocked_pass_all(:,2)==12); % select cycles that passed containing 12 upd
blocked_pass_12=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_12.dat',blocked_pass_12, ' ')
indcs=(blocked_fail_all(:,2)==12); % select cycles that failed containing 12 upd
blocked_fail_12=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_12.dat',blocked_fail_12, ' ')

indcs=(blocked_pass_all(:,2)==13); % select cycles that passed containing 13 upd
blocked_pass_13=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_13.dat',blocked_pass_13, ' ')
indcs=(blocked_fail_all(:,2)==13); % select cycles that failed containing 13 upd
blocked_fail_13=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_13.dat',blocked_fail_13, ' ')

indcs=(blocked_pass_all(:,2)==14); % select cycles that passed containing 14 upd
blocked_pass_14=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_14.dat',blocked_pass_14, ' ')
indcs=(blocked_fail_all(:,2)==14); % select cycles that failed containing 14 upd
blocked_fail_14=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_14.dat',blocked_fail_14, ' ')

indcs=(blocked_pass_all(:,2)==15); % select cycles that passed containing 15 upd
blocked_pass_15=blocked_pass_all(indcs,:);
dlmwrite('blocked_pass_15.dat',blocked_pass_15, ' ')
indcs=(blocked_fail_all(:,2)==15); % select cycles that failed containing 15 upd
blocked_fail_15=blocked_fail_all(indcs,:);
dlmwrite('blocked_fail_15.dat',blocked_fail_15, ' ')


indcs=(data_lapack(:,5)==0); % select cycles that passed
lapack_pass_all=data_lapack(indcs,:);
dlmwrite('lapack_pass_all.dat',lapack_pass_all, ' ')

indcs=(data_lapack(:,5)!=0); % select cycles that failed
lapack_fail_all=data_lapack(indcs,:);
dlmwrite('lapack_fail_all.dat',lapack_fail_all, ' ')

indcs=(lapack_pass_all(:,2)==1); % select cycles that passed containing 1 upd
lapack_pass_1=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_1.dat',lapack_pass_1, ' ')
indcs=(lapack_fail_all(:,2)==1); % select cycles that failed containing 1 upd
lapack_fail_1=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_1.dat',lapack_fail_1, ' ')

indcs=(lapack_pass_all(:,2)==2); % select cycles that passed containing 2 upd
lapack_pass_2=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_2.dat',lapack_pass_2, ' ')
indcs=(lapack_fail_all(:,2)==2); % select cycles that failed containing 2 upd
lapack_fail_2=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_2.dat',lapack_fail_2, ' ')

indcs=(lapack_pass_all(:,2)==3); % select cycles that passed containing 3 upd
lapack_pass_3=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_3.dat',lapack_pass_3, ' ')
indcs=(lapack_fail_all(:,2)==3); % select cycles that failed containing 3 upd
lapack_fail_3=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_3.dat',lapack_fail_3, ' ')

indcs=(lapack_pass_all(:,2)==4); % select cycles that passed containing 4 upd
lapack_pass_4=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_4.dat',lapack_pass_4, ' ')
indcs=(lapack_fail_all(:,2)==4); % select cycles that failed containing 4 upd
lapack_fail_4=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_4.dat',lapack_fail_4, ' ')

indcs=(lapack_pass_all(:,2)==5); % select cycles that passed containing 5 upd
lapack_pass_5=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_5.dat',lapack_pass_5, ' ')
indcs=(lapack_fail_all(:,2)==5); % select cycles that failed containing 5 upd
lapack_fail_5=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_5.dat',lapack_fail_5, ' ')

indcs=(lapack_pass_all(:,2)==6); % select cycles that passed containing 6 upd
lapack_pass_6=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_6.dat',lapack_pass_6, ' ')
indcs=(lapack_fail_all(:,2)==6); % select cycles that failed containing 6 upd
lapack_fail_6=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_6.dat',lapack_fail_6, ' ')

indcs=(lapack_pass_all(:,2)==7); % select cycles that passed containing 7 upd
lapack_pass_7=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_7.dat',lapack_pass_7, ' ')
indcs=(lapack_fail_all(:,2)==7); % select cycles that failed containing 7 upd
lapack_fail_7=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_7.dat',lapack_fail_7, ' ')

indcs=(lapack_pass_all(:,2)==8); % select cycles that passed containing 8 upd
lapack_pass_8=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_8.dat',lapack_pass_8, ' ')
indcs=(lapack_fail_all(:,2)==8); % select cycles that failed containing 8 upd
lapack_fail_8=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_8.dat',lapack_fail_8, ' ')

indcs=(lapack_pass_all(:,2)==9); % select cycles that passed containing 9 upd
lapack_pass_9=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_9.dat',lapack_pass_9, ' ')
indcs=(lapack_fail_all(:,2)==9); % select cycles that failed containing 9 upd
lapack_fail_9=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_9.dat',lapack_fail_9, ' ')

indcs=(lapack_pass_all(:,2)==10); % select cycles that passed containing 10 upd
lapack_pass_10=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_10.dat',lapack_pass_10, ' ')
indcs=(lapack_fail_all(:,2)==10); % select cycles that failed containing 10 upd
lapack_fail_10=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_10.dat',lapack_fail_10, ' ')

indcs=(lapack_pass_all(:,2)==11); % select cycles that passed containing 11 upd
lapack_pass_11=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_11.dat',lapack_pass_11, ' ')
indcs=(lapack_fail_all(:,2)==11); % select cycles that failed containing 11 upd
lapack_fail_11=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_11.dat',lapack_fail_11, ' ')

indcs=(lapack_pass_all(:,2)==12); % select cycles that passed containing 12 upd
lapack_pass_12=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_12.dat',lapack_pass_12, ' ')
indcs=(lapack_fail_all(:,2)==12); % select cycles that failed containing 12 upd
lapack_fail_12=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_12.dat',lapack_fail_12, ' ')

indcs=(lapack_pass_all(:,2)==13); % select cycles that passed containing 13 upd
lapack_pass_13=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_13.dat',lapack_pass_13, ' ')
indcs=(lapack_fail_all(:,2)==13); % select cycles that failed containing 13 upd
lapack_fail_13=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_13.dat',lapack_fail_13, ' ')

indcs=(lapack_pass_all(:,2)==14); % select cycles that passed containing 14 upd
lapack_pass_14=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_14.dat',lapack_pass_14, ' ')
indcs=(lapack_fail_all(:,2)==14); % select cycles that failed containing 14 upd
lapack_fail_14=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_14.dat',lapack_fail_14, ' ')

indcs=(lapack_pass_all(:,2)==15); % select cycles that passed containing 15 upd
lapack_pass_15=lapack_pass_all(indcs,:);
dlmwrite('lapack_pass_15.dat',lapack_pass_15, ' ')
indcs=(lapack_fail_all(:,2)==15); % select cycles that failed containing 15 upd
lapack_fail_15=lapack_fail_all(indcs,:);
dlmwrite('lapack_fail_15.dat',lapack_fail_15, ' ')


indcs=(data_wb2(:,5)==0); % select cycles that passed
wb2_pass_all=data_wb2(indcs,:);
dlmwrite('wb2_pass_all.dat',wb2_pass_all, ' ')

indcs=(data_wb2(:,5)!=0); % select cycles that failed
wb2_fail_all=data_wb2(indcs,:);
dlmwrite('wb2_fail_all.dat',wb2_fail_all, ' ')


indcs=(data_wb3(:,5)==0); % select cycles that passed
wb3_pass_all=data_wb3(indcs,:);
dlmwrite('wb3_pass_all.dat',wb3_pass_all, ' ')

indcs=(data_wb3(:,5)!=0); % select cycles that failed
wb3_fail_all=data_wb3(indcs,:);
dlmwrite('wb3_fail_all.dat',wb3_fail_all, ' ')


n_all_cycles=size(data_anthony)(1);
n_1_cycles=size(anthony_pass_1)(1)+size(anthony_fail_1)(1);
n_2_cycles=size(anthony_pass_2)(1)+size(anthony_fail_2)(1);
n_3_cycles=size(anthony_pass_3)(1)+size(anthony_fail_3)(1);
n_6_cycles=size(anthony_pass_6)(1)+size(anthony_fail_6)(1);

fail_rate_all_anthony=sum(anthony_fail_all(:,5))/n_all_cycles;
fail_rate_all_naive=sum(naive_fail_all(:,5))/n_all_cycles;
fail_rate_all_later=sum(later_fail_all(:,5))/n_all_cycles;
fail_rate_all_split=sum(split_fail_all(:,5))/n_all_cycles;
fail_rate_all_blocked=sum(blocked_fail_all(:,5))/n_all_cycles;
fail_rate_all_lapack=sum(lapack_fail_all(:,5))/n_all_cycles;

fail_rate_1_anthony=sum(anthony_fail_1(:,5))/n_1_cycles;
fail_rate_1_naive=sum(naive_fail_1(:,5))/n_1_cycles;
fail_rate_1_later=sum(later_fail_1(:,5))/n_1_cycles;
fail_rate_1_split=sum(split_fail_1(:,5))/n_1_cycles;
fail_rate_1_blocked=sum(blocked_fail_1(:,5))/n_1_cycles;
fail_rate_1_lapack=sum(lapack_fail_1(:,5))/n_1_cycles;

fail_rate_2_anthony=sum(anthony_fail_2(:,5))/n_2_cycles;
fail_rate_2_naive=sum(naive_fail_2(:,5))/n_2_cycles;
fail_rate_2_later=sum(later_fail_2(:,5))/n_2_cycles;
fail_rate_2_split=sum(split_fail_2(:,5))/n_2_cycles;
fail_rate_2_blocked=sum(blocked_fail_2(:,5))/n_2_cycles;
fail_rate_2_lapack=sum(lapack_fail_2(:,5))/n_2_cycles;
fail_rate_wb2=sum(data_wb2(:,5))/n_2_cycles;

fail_rate_3_anthony=sum(anthony_fail_3(:,5))/n_3_cycles;
fail_rate_3_naive=sum(naive_fail_3(:,5))/n_3_cycles;
fail_rate_3_later=sum(later_fail_3(:,5))/n_3_cycles;
fail_rate_3_split=sum(split_fail_3(:,5))/n_3_cycles;
fail_rate_3_blocked=sum(blocked_fail_3(:,5))/n_3_cycles;
fail_rate_3_lapack=sum(lapack_fail_3(:,5))/n_3_cycles;
fail_rate_wb3=sum(data_wb3(:,5))/n_3_cycles;

fail_rate_6_anthony=sum(anthony_fail_6(:,5))/n_6_cycles;
fail_rate_6_naive=sum(naive_fail_6(:,5))/n_6_cycles;
fail_rate_6_later=sum(later_fail_6(:,5))/n_6_cycles;
fail_rate_6_split=sum(split_fail_6(:,5))/n_6_cycles;
fail_rate_6_lapack=sum(lapack_fail_6(:,5))/n_6_cycles;
fail_rate_6_blocked=sum(blocked_fail_6(:,5))/n_6_cycles;

printf("\n");
printf ("Fail rates for all (N=%d) cycles\n", n_all_cycles);
printf ("-------------------------------------------------------------------------------------------\n");
printf ("Anthony:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_all_anthony*100, fail_rate_all_anthony, size(anthony_pass_all)(1), size(anthony_fail_all)(1), size(anthony_pass_all)(1)+size(anthony_fail_all)(1));
printf ("Naive:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_all_naive*100, fail_rate_all_naive, size(naive_pass_all)(1), size(naive_fail_all)(1), size(naive_pass_all)(1)+size(naive_fail_all)(1));
printf ("Later:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_all_later*100, fail_rate_all_later, size(later_pass_all)(1), size(later_fail_all)(1), size(later_pass_all)(1)+size(later_fail_all)(1));
printf ("Splitting:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_all_split*100, fail_rate_all_split, size(split_pass_all)(1), size(split_fail_all)(1), size(split_pass_all)(1)+size(split_fail_all)(1));
printf ("Blocked:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_all_blocked*100, fail_rate_all_blocked, size(blocked_pass_all)(1), size(blocked_fail_all)(1), size(blocked_pass_all)(1)+size(blocked_fail_all)(1));
printf ("Lapack:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_all_lapack*100, fail_rate_all_lapack, size(lapack_pass_all)(1), size(lapack_fail_all)(1), size(lapack_pass_all)(1)+size(lapack_fail_all)(1));

printf("\n");
printf ("Fail rates for cycles containing 1 update (N=%d) (solely due to numerical noise)\n", n_1_cycles);
printf ("-------------------------------------------------------------------------------------------\n");
printf ("Anthony:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_1_anthony*100, fail_rate_1_anthony, size(anthony_pass_1)(1), size(anthony_fail_1)(1), n_1_cycles);
printf ("Naive:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_1_naive*100, fail_rate_1_naive, size(naive_pass_1)(1), size(naive_fail_1)(1), size(naive_pass_1)(1)+size(naive_fail_1)(1));
printf ("Later:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_1_later*100, fail_rate_1_later, size(later_pass_1)(1), size(later_fail_1)(1), size(later_pass_1)(1)+size(later_fail_1)(1));
printf ("Splitting:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_1_split*100, fail_rate_1_split, size(split_pass_1)(1), size(split_fail_1)(1), size(split_pass_1)(1)+size(split_fail_1)(1));
printf ("Blocked:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_1_blocked*100, fail_rate_1_blocked, size(blocked_pass_1)(1), size(blocked_fail_1)(1), size(blocked_pass_1)(1)+size(blocked_fail_1)(1));
printf ("Lapack:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_1_lapack*100, fail_rate_1_lapack, size(lapack_pass_1)(1), size(lapack_fail_1)(1), size(lapack_pass_1)(1)+size(lapack_fail_1)(1));

printf("\n");
printf ("Fail rates for cycles containing 2 updates (N=%d) (compare blocked w/ WB2)\n", n_2_cycles);
printf ("-------------------------------------------------------------------------------------------\n");
printf ("Anthony:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_2_anthony*100, fail_rate_2_anthony, size(anthony_pass_2)(1), size(anthony_fail_2)(1), n_2_cycles);
printf ("Naive:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_2_naive*100, fail_rate_2_naive, size(naive_pass_2)(1), size(naive_fail_2)(1), size(naive_pass_2)(1)+size(naive_fail_2)(1));
printf ("Later:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_2_later*100, fail_rate_2_later, size(later_pass_2)(1), size(later_fail_2)(1), size(later_pass_2)(1)+size(later_fail_2)(1));
printf ("Splitting:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_2_split*100, fail_rate_2_split, size(split_pass_2)(1), size(split_fail_2)(1), size(split_pass_2)(1)+size(split_fail_2)(1));
printf ("Blocked:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_2_blocked*100, fail_rate_2_blocked, size(blocked_pass_2)(1), size(blocked_fail_2)(1), size(blocked_pass_2)(1)+size(blocked_fail_2)(1));
printf ("Lapack:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_2_lapack*100, fail_rate_2_lapack, size(lapack_pass_2)(1), size(lapack_fail_2)(1), size(lapack_pass_2)(1)+size(lapack_fail_2)(1));
printf ("Woodbury 2:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_wb2*100, fail_rate_wb2, size(wb2_pass_all)(1), size(wb2_fail_all)(1), size(wb2_pass_all)(1)+size(wb2_fail_all)(1));

printf("\n");
printf ("Fail rates for cycles containing 3 updates (N=%d) (compare blocked w/ WB3)\n", n_3_cycles);
printf ("-------------------------------------------------------------------------------------------\n");
printf ("Anthony:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_3_anthony*100, fail_rate_3_anthony, size(anthony_pass_3)(1), size(anthony_fail_3)(1), n_3_cycles);
printf ("Naive:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_3_naive*100, fail_rate_3_naive, size(naive_pass_3)(1), size(naive_fail_3)(1), size(naive_pass_3)(1)+size(naive_fail_3)(1));
printf ("Later:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_3_later*100, fail_rate_3_later, size(later_pass_3)(1), size(later_fail_3)(1), size(later_pass_3)(1)+size(later_fail_3)(1));
printf ("Splitting:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_3_split*100, fail_rate_3_split, size(split_pass_3)(1), size(split_fail_3)(1), size(split_pass_3)(1)+size(split_fail_3)(1));
printf ("Blocked:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_3_blocked*100, fail_rate_3_blocked, size(blocked_pass_3)(1), size(blocked_fail_3)(1), size(blocked_pass_3)(1)+size(blocked_fail_3)(1));
printf ("Lapack:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_3_lapack*100, fail_rate_3_lapack, size(lapack_pass_3)(1), size(lapack_fail_3)(1), size(lapack_pass_3)(1)+size(lapack_fail_3)(1));
printf ("Woodbury 3:\t%f (= %f x N cycles; %d pass + %d fail = %d tot.)\n", fail_rate_wb3*100, fail_rate_wb3, size(wb3_pass_all)(1), size(wb3_fail_all)(1), size(wb3_pass_all)(1)+size(wb3_fail_all)(1));

printf("\n");
printf ("Fail rates for cycles containing 6 updates (N=%d) (blocked vs splitting in multiples of 3)\n", n_6_cycles);
printf ("-------------------------------------------------------------------------------------------\n");
printf ("Anthony:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_6_anthony*100, fail_rate_6_anthony, size(anthony_pass_6)(1), size(anthony_fail_6)(1), n_6_cycles);
printf ("Naive:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_6_naive*100, fail_rate_6_naive, size(naive_pass_6)(1), size(naive_fail_6)(1), size(naive_pass_6)(1)+size(naive_fail_6)(1));
printf ("Later:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_6_later*100, fail_rate_6_later, size(later_pass_6)(1), size(later_fail_6)(1), size(later_pass_6)(1)+size(later_fail_6)(1));
printf ("Splitting:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_6_split*100, fail_rate_6_split, size(split_pass_6)(1), size(split_fail_6)(1), size(split_pass_6)(1)+size(split_fail_6)(1));
printf ("Blocked:\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_6_blocked*100, fail_rate_6_blocked, size(blocked_pass_6)(1), size(blocked_fail_6)(1), size(blocked_pass_6)(1)+size(blocked_fail_6)(1));
printf ("Lapack:\t\t%f (= %f x N cycles; %d pass + %d fail = %d tot)\n", fail_rate_6_lapack*100, fail_rate_6_lapack, size(lapack_pass_6)(1), size(lapack_fail_6)(1), size(lapack_pass_6)(1)+size(lapack_fail_6)(1));
