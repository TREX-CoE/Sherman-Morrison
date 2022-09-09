void print_m(const double* mat, uint16_t m, uint16_t n, uint16_t ldm, char* name)
{
  printf("%s = \n", name);
  for (uint16_t i = 0; i < m; ++i)
  {
    for (uint16_t j = 0; j < n; ++j)
    {
      printf("%11.5f ", mat[i * ldm + j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_m_t(const double* mat, uint16_t m, uint16_t n, uint16_t ldm, char* name)
{
  printf("%s = \n", name);
  for (uint16_t i = 0; i < m; ++i)
  {
    for (uint16_t j = 0; j < n; ++j)
    {
      printf("%11.5f ", mat[j * ldm + i]);
    }
    printf("\n");
  }
  printf("\n");
}

void transpose(double* a, uint16_t lda, double *b, uint16_t ldb, uint16_t m, uint16_t n)
{
  for(uint16_t i = 0; i < m; i++)
  {
    for( uint16_t j = 0; j < n; j++)
    {
      b[j * ldb + i] = a[i * lda + j];
    }
  }
}
