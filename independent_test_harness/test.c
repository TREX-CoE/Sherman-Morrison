#include <stdint.h>
#include "meuk.h"
#include "cycles.h"
#include "debug.h"

#define DATASET "dataset"
#define REPETITIONS 1

uint64_t n_splits;
uint64_t block_fail;
uint64_t recursive_calls;

int main(int argc, char **argv) {
  assert(argc == 2);
  char *version = argv[1];

#ifdef HAVE_CUBLAS_OFFLOAD
  cublasHandle_t handle = init_cublas();
  cusolverDnHandle_t s_handle = init_cusolver();
#endif

  // SETUP STORAGE AND DATA ACCESS
  hid_t  file_id = H5Fopen(DATASET, H5F_ACC_RDONLY, H5P_DEFAULT);

  printf("#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  printf("#1\t2\t3\t4\t\t5\t6\t\t7\t\t8\t\t9\t\t10\t\t11\t\t12\t\t13\t\t14\n");
  printf("#CYCLE\tUPDS\tERR_IN\tERR_BREAK\tERR_OUT\tSPLITS\t\tBLK_FAILS\tMAX\t\tFROB\t\tCOND\t\tCPU_CYC\t\tCPU_CYC/UPD\tCUMUL\t\tREC\n");
  printf("#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  // FOR EACH UPDATE CYCLE DO:
  // for (uint32_t cycles_index = 0; cycles_index < n_cycles; cycles_index++) {
  for (uint32_t cycles_index = 0; cycles_index < 1; cycles_index++) {

    // SETUP TEST PARAMETERS
    const double breakdown = 0.001; // default = 0.001. 1e-9 might be too small
    const double tolerance = 0.001; // default = 0.001
    double cumulative = 0;

    // 1. READ DATA FROM DATASET
    uint32_t cycle = cycles[cycles_index];
    uint64_t Dim = get_dim(cycle, file_id);
    uint64_t Lds = Dim;
    uint64_t N_updates = get_nupdates(cycle, file_id);
    uint64_t* Updates_index = get_upd_idcs(cycle, file_id, N_updates);
    double* Updates = get_upds(cycle, file_id, N_updates, Lds);
    double* Slater = get_slater(cycle, file_id, Dim, Lds);
    double* Slater_invT = get_slater_inv(cycle, file_id, Dim, Lds);
    double *Slater_invT_copy = calloc(1, sizeof *Slater_invT_copy * Dim * Lds);
    double determinant = get_determinant(cycle, file_id), determinant_copy;

    // Compute transpose of S. ST: 24 x 21
    double *SlaterT = calloc(1, sizeof *SlaterT * Dim * Lds);
    transpose(Slater, Lds, SlaterT, Dim, Dim, Lds);

    // Convert repl. upds into additive upds.
    convert(Updates, N_updates, Updates_index, SlaterT, Dim, Lds);

    // 2. CHECK ERROR ON THE INPUT DATA AND RECORD RESULT: ERR_INPUT
    uint32_t err_inp = check_error(Lds, Dim, Slater_invT, SlaterT, tolerance);

    // Update Slater matrix
    update(SlaterT, Updates, Updates_index, N_updates, Dim, Lds);

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
      memcpy(Slater_invT_copy, Slater_invT, Lds * Dim * sizeof(double));
      determinant_copy = determinant;

      // ### CHOOSE A KERNEL:
      if (version[0] == 'n') { // Naive
        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_sherman_morrison(Lds, Dim, N_updates, Updates,
                                           Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
      } else if (version[0] == 'l') { // Later

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_sherman_morrison_later(Lds, Dim, N_updates, Updates,
                                                 Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
      } else if (version[0] == '2') { // by twos

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_woodbury_2(Lds, Dim, Updates, Updates_index,
                                     breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);

      } else if (version[0] == '3') { // by threes

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_woodbury_3(Lds, Dim, Updates, Updates_index,
                                     breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);

      } else if (version[0] == 'k') { // Woodbury K

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_woodbury_k(Lds, Dim, N_updates, Updates,
                                     Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);

#ifdef HAVE_CUBLAS_OFFLOAD
      } else if (version[0] == 'c') { // Woodbury K cuBLAS

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_woodbury_k_cublas_offload(handle, s_handle, Lds, Dim, N_updates, Updates,
                                                    Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
#endif
      } else if (version[0] == 's') { // Splitting

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_sherman_morrison_splitting(Lds, Dim, N_updates, Updates,
                                                     Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
      } else if (version[0] == 'b') { // Blocked

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = qmckl_sherman_morrison_smw32s(Lds, Dim, N_updates, Updates,
                                                  Updates_index, breakdown, Slater_invT_copy, &determinant);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);
      } else if (version[0] == 'm') { // LAPACK/MKL

        // Only send upper Dim x Dim part of matrix to lapack
        double *tmp = malloc(sizeof *tmp * Dim * Dim);
        memcpy(tmp, SlaterT, sizeof *SlaterT * Dim * Dim);

        // 1. FETCH START TIME
        uint64_t before = rdtsc();

        // 2. EXECUTE KERNEL AND REMEMBER EXIT STATUS
        err_break = inverse(tmp, Dim, Dim);

        // 3. FETCH FINISH TIME
        uint64_t after = rdtsc();

        // Copy elements of inverse back, adding 0-padding in "correct" place
        copy(Slater_invT_copy, Lds, tmp, Dim);
        free(tmp);

        // 4. ADD TIME DIFFERENCE TO TIME CUMMULATOR
        accumulator += (double)(after - before);

      } else { // Exit
        printf("Version '%c' not implemented.\n", version[0]);
        return 1;
      }

    } // END OF REPETITIONS LOOP

    // 4. COPY RESULT BACK TO ORIGINAL
    memcpy(Slater_invT, Slater_invT_copy, Lds * Dim * sizeof(double));
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

    // Compute Slater x Slater_inv
    double* SSi = malloc(sizeof *SSi * Dim * Dim);
    double alpha = 1.0, beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                Dim, Dim, Dim,
                alpha, SlaterT, Dim, Slater_invT, Lds,
                beta, SSi, Dim);

    // Compute residual matrix S * Sinv - Id
    double* Res = calloc(1, sizeof *Res * Dim * Dim);
    residual(SSi, Res, Dim);
    const double max = max_norm(Res, Dim, Dim);
    free(SSi);

    // 7. CHECK ERRROR ON THE UPDATED DATA AND RECORD THE RESULT: ERR_OUT
    uint32_t err_out = check_error(Lds, Dim, Slater_invT, SlaterT, tolerance);
    // int32_t err_out = check_error_better(max, tolerance);

    // if (err_out == 1) printf("cycle index %d: cycle %d with %lu upds failed!\n", cycles_index, cycle, N_updates);

    // 8. COMPUTE CONDITION NUMBER
    const double condnr = condition_number(Slater, Slater_invT, Lds, Dim);
    const double frob = frobenius_norm(Res, Dim, Dim);
    free(Res);

    // 10. WRITE RESULTS TO FILE: CYCLE#, #UPDS, ERR_INP, ERR_BREAK, #SPLITS, ERR_OUT, COND, #CLCK_TCKS
    printf("%u\t%lu\t%u\t%u\t\t%u\t%lu\t\t%lu\t\t%e\t%e\t%e\t%e\t%e\t%e\t%lu\n", cycle, N_updates, err_inp, err_break, err_out, n_splits, block_fail, max, frob, condnr, accumulator, cycles_per_update, cumulative, recursive_calls);

    free(Updates_index);
    free(Updates);
    free(SlaterT);
    free(Slater_invT);
    free(Slater_invT_copy);

  } // END OF CYCLE LOOP
  printf("#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  printf("#1\t2\t3\t4\t\t5\t6\t\t7\t\t8\t\t9\t\t10\t\t11\t\t12\t\t13\t\t14\n");
  printf("#CYCLE\tUPDS\tERR_IN\tERR_BREAK\tERR_OUT\tSPLITS\t\tBLK_FAILS\tMAX\t\tFROB\t\tCOND\t\tCPU_CYC\t\tCPU_CYC/UPD\tCUMUL\t\tREC\n");
  printf("#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");

  (void) H5Fclose(file_id);

#ifdef HAVE_CUBLAS_OFFLOAD
  cublasDestroy(handle);
  cusolverDnDestroy(s_handle);
#endif
}
