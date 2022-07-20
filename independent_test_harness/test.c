#include "meuk.h"
#include "cycles.h"
#include <stdint.h>

#define DATASET "dataset_329d_zeropadded_cm.hdf5"
// #define DATASET "dataset_15784d_zeropadded_cm.hdf5"
#define REPETITIONS 100000

uint64_t n_splits;
uint64_t block_fail;
uint64_t recursive_calls;

int main(int argc, char **argv) {
  assert(argc == 2);
  char *version = argv[1];

  // SETUP STORAGE AND DATA ACCESS
  hid_t  file_id, dataset_id;
  herr_t status;
  file_id = H5Fopen(DATASET, H5F_ACC_RDONLY, H5P_DEFAULT);
  char nupds_key[32];
  char upd_idx_key[32];
  char upds_key[32];
  char slater_key[32];
  char slater_inv_key[32];
  char det_key[32];
  const uint64_t Dim = 21;
  const uint64_t LDS = 24;
  uint64_t N_updates;
  double Slater[LDS * Dim ], SlaterT[LDS * Dim];
  double Slater_invT[LDS * Dim], Slater_invT_copy[LDS * Dim];
  double determinant, determinant_copy;

  // SETUP TEST PARAMETERS
  const double breakdown = 0.001; // default = 0.001. 1e-9 might be too small
  const double tolerance = 0.001; // default = 0.001
  double cumulative = 0;

printf("#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
printf("#1\t2\t3\t4\t\t5\t6\t\t7\t\t8\t\t9\t\t10\t\t11\t\t12\t\t13\t\t14\n");
printf("#CYCLE\tUPDS\tERR_IN\tERR_BREAK\tERR_OUT\tSPLITS\t\tBLK_FAILS\tMAX\t\tFROB\t\tCOND\t\tCPU_CYC\t\tCPU_CYC/UPD\tCUMUL\t\tREC\n");
printf("#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  // FOR EACH UPDATE CYCLE DO:
  // for (uint32_t cycles_index = 0; cycles_index < n_cycles; cycles_index++) {
  for (uint32_t cycles_index = 40; cycles_index < 41; cycles_index++) {  
    // 1. READ DATA FROM DATASET
    uint32_t cycle = cycles[cycles_index];
    sprintf(nupds_key, "/cycle_%d/nupdates", cycle);
    sprintf(upd_idx_key, "/cycle_%d/col_update_index", cycle);
    sprintf(upds_key, "/cycle_%d/updates", cycle);
    sprintf(slater_key, "/cycle_%d/slater_matrix", cycle);
    sprintf(slater_inv_key, "/cycle_%d/slater_inverse_t", cycle);
    sprintf(det_key, "/cycle_%d/determinant", cycle);
    read_uint(file_id, nupds_key, &N_updates);
    uint64_t *Updates_index = malloc(N_updates * sizeof(uint64_t));
    double *Updates = malloc(LDS * N_updates * sizeof(double));
    read_uint(file_id, upd_idx_key, Updates_index);
    read_double(file_id, upds_key, Updates);
    read_double(file_id, slater_key, Slater);
    read_double(file_id, slater_inv_key, Slater_invT);
    read_double(file_id, det_key, &determinant);

    // Compute transpose of S. ST: 24 x 21
    for (int i = 0; i < LDS; i++) {
      for (int j = 0; j < Dim; j++) {
        SlaterT[i * Dim + j] = Slater[j * LDS + i];
      }
    }

    // Convert repl. upds into additive upds.
    for (int i = 0; i < N_updates; i++) {
      int col = Updates_index[i] - 1;
      for (int j = 0; j < LDS; j++) {
        Updates[i * LDS + j] -= SlaterT[col + j * Dim];
      }
    }

    // 2. CHECK ERROR ON THE INPUT DATA AND RECORD RESULT: ERR_INPUT
    uint32_t err_inp = check_error(LDS, Dim, Slater_invT, SlaterT, tolerance);

    // Update Slater matrix
    for (int i = 0; i < N_updates; i++) {
      int col = Updates_index[i] - 1;
      for (int j = 0; j < Dim; j++) {
        SlaterT[col + j * Dim] += Updates[i * LDS + j];
      }
    } // A this point SlaterT, Updates & the updated SlaterT are correct. Checked in GDB

    int32_t err_break;

    // 3. SET TIME- AND SPLIT ACCUMULATOR TO ZERO
    double accumulator = 0;
    double cycles_per_update = 0;
    n_splits = 0;
    block_fail = 0;
    recursive_calls = 0;

    // ## FOR A SET NUMBER OF REPETITIONS DO:
    for (int rep = 0; rep < REPETITIONS; rep++) {

      // 1. MAKE A FRESH COPY OF THE SLATER INVERSE AND DETERMINANT AND USE THE COPY
      memcpy(Slater_invT_copy, Slater_invT, LDS * Dim * sizeof(double));
      determinant_copy = determinant;

      // ### CHOOSE A KERNEL:
      if (version[0] == 'a') { // Anthony
        const double *Upds;
        const uint64_t *Ui;
        double determinant_previous;

        err_break = 0;

        for (int i = 0; i < LDS * Dim; i++) Slater_invT_copy[i] *= determinant_copy; // Multiply inv(Slater-mat) by det(Slater-mat) to get adj(Slater_mat)

        for (int i = 0; i < N_updates; i++) {
          Upds = &Updates[i * LDS];
          Ui = &Updates_index[i];
          determinant_previous = determinant_copy;

          // 1. FETCH START TIME
          uint64_t before = rdtsc();

          // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
          detupd(Dim, LDS, Upds, Ui, Slater_invT_copy, &determinant_copy);

          // 3. FETCH FINISH TIME
          uint64_t after = rdtsc();

          // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
          accumulator += (double)(after - before);

          // 5. STOP APPLYING UPDATES IF BREAKDOWN DETECTED
          double lambda = determinant_copy / determinant_previous; // should be id. to lambda in detupd
          if (fabs(lambda) < breakdown) {
            err_break = 1;
            break;
          }
        }
        
        if (err_break == 1) { // Divide adj(Slater-mat) by OLD det(Slater-mat) to get inv(Slater_mat) again
          for (int i = 0; i < LDS * Dim; i++) Slater_invT_copy[i] /= determinant_previous;
        } else { // Divide adj(Slater-mat) by NEW det(Slater-mat) to get inv(Slater_mat) again
          for (int i = 0; i < LDS * Dim; i++) Slater_invT_copy[i] /= determinant_copy;
        }
      } else if (version[0] == 'n') { // Naive

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_sherman_morrison(LDS, Dim, N_updates, Updates,
              Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
      } else if (version[0] == 'l') { // Later

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_sherman_morrison_later(LDS, Dim, N_updates, Updates,
              Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
      } else if (version[0] == '2') { // by twos

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_woodbury_2(LDS, Dim, Updates, Updates_index,
              breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);

      } else if (version[0] == '3') { // by threes

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_woodbury_3(LDS, Dim, Updates, Updates_index,
              breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);

      } else if (version[0] == 'k') { // Woodbury K

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_woodbury_k(LDS, Dim, N_updates, Updates,
              Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);

      } else if (version[0] == 's') { // Splitting

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_sherman_morrison_splitting(LDS, Dim, N_updates, Updates,
              Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
      } else if (version[0] == 'b') { // Blocked

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_sherman_morrison_smw32s(LDS, Dim, N_updates, Updates,
              Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
      } else if (version[0] == 'm') { // LAPACK/MKL

        // Only send upper Dim x Dim part of matrix to lapack
        double tmp[Dim*Dim];
        memcpy(tmp, SlaterT, Dim*Dim*sizeof(double));

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = inverse(tmp, Dim, Dim);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // Copy elements of inverse back, adding 0-padding in "correct" place
        for (uint32_t i = 0; i < Dim; i++) {
          for (uint32_t j = 0; j < LDS; j++) {
            if (j < Dim) Slater_invT_copy[i * LDS + j] = tmp[i * Dim + j];
            else Slater_invT_copy[i * LDS + j] = 0.0;
          }
        }

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);

      } else {                        // Exit
        printf("Version '%c' not implemented.\n", version[0]);
        return 1;
      }
    } // END OF REPETITIONS LOOP

    // 4. COPY RESULT BACK TO ORIGINAL 
    memcpy(Slater_invT, Slater_invT_copy, LDS * Dim * sizeof(double));
    determinant = determinant_copy;
    // At this point Slater_invT contains the correct inverse matrix

    // 5. DIVIDE CYCLE- AND SPLIT-ACCUMULATOR BY NUMBER OF REPETITIONS AND RECORD
    //    DIVIDE CYCLE-ACCUMULATOR BY NUMBER OF UPDATES AND RECORD
    accumulator /= REPETITIONS;
    cycles_per_update = accumulator / N_updates;
    n_splits /= REPETITIONS;
    block_fail /= REPETITIONS;
    recursive_calls /= REPETITIONS;

    // 6. ADD THE AVERAGED TIME PER CYCLE OF ACCUMULATER TO
    // CUMULATIVE RESULT FOR THE ENTIRE DATASET
    cumulative += accumulator;

    double SSi[Dim * Dim];
    matmul(SlaterT, Slater_invT, SSi, LDS, Dim);
    double Res[Dim * Dim];
    residual(SSi, Res, Dim);
    const double max = max_norm(Res, Dim, Dim);

    // 7. CHECK ERRROR ON THE UPDATED DATA AND RECORD THE RESULT: ERR_OUT
    uint32_t err_out = check_error(LDS, Dim, Slater_invT, SlaterT, tolerance);
    // int32_t err_out = check_error_better(max, tolerance);

    // if (err_out == 1) printf("cycle index %d: cycle %d with %lu upds failed!\n", cycles_index, cycle, N_updates);

    // 8. COMPUTE CONDITION NUMBER
    const double condnr = condition_number(Slater, Slater_invT, LDS, Dim);
    const double frob = frobenius_norm(Res, Dim, Dim);


    // 10. WRITE RESULTS TO FILE: CYCLE#, #UPDS, ERR_INP, ERR_BREAK, #SPLITS, ERR_OUT, COND, #CLCK_TCKS
        printf("%u\t%lu\t%u\t%u\t\t%u\t%lu\t\t%lu\t\t%e\t%e\t%e\t%e\t%e\t%e\t%lu\n", cycle, N_updates, err_inp, err_break, err_out, n_splits, block_fail, max, frob, condnr, accumulator, cycles_per_update, cumulative, recursive_calls);
        
    free(Updates_index);
    free(Updates);

  } // END OF CYCLE LOOP
  printf("#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  printf("#1\t2\t3\t4\t\t5\t6\t\t7\t\t8\t\t9\t\t10\t\t11\t\t12\t\t13\t\t14\n");
  printf("#CYCLE\tUPDS\tERR_IN\tERR_BREAK\tERR_OUT\tSPLITS\t\tBLK_FAILS\tMAX\t\tFROB\t\tCOND\t\tCPU_CYC\t\tCPU_CYC/UPD\tCUMUL\t\tREC\n");
  printf("#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");

  (void) H5Fclose(file_id);
}
