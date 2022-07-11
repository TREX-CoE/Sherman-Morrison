Test method
===========

# > for each update cycle do:
    (# of updates changes -> update indices & size of update-matrix changes)

    1. read data from dataset
    2. check error on the input data and record result: ERR_INPUT
    3. set cycle- and split accumulator to zero
    
 ## >> for a set number of repetitions do:

        1. take a fresh copy (memcpy) of the slater inverse and use it in chosen kernel

 ### >>> for the chosen kernel do:

            1. fetch start cycles
            2. execute kernel and remember exit status: ERR_BREAK
               (number of splits is recorded in global variable) 
            3. fetch finish cycles
            4. add cycle difference to time acummulator

## > continue: for each update cycle do

    4. copy the updated slater-inverse-copy back to original
    5a. divide cycle- and split-accumulator by number of repetitions
    5b. divide cycle-accumulator by number of updates
    6. add the averaged time/update-cycle of accumulater to cummulative-
        result for the entire dataset
    7. update the slater matrix
    8. check the error on the updated data and record the result: ERR_OUT
    9. write results to stdout: cycle#, #upds, err_inp, err_break, 
        #splits, err_out, #clck_tcks, #clck_tcks/upd, cumulative cycles
