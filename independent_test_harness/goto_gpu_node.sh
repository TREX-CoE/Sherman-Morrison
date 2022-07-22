salloc --nodes=1 --account=prcoe10 -p booster --gres=gpu:1
wait
srun --pty /bin/bash

