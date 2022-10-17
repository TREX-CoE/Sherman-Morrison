#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

void print_dm(const double* mat, uint16_t m, uint16_t n, uint16_t ldm, char* name)
{
  printf("%s = \n", name);
  for (uint16_t i = 0; i < m; ++i)
  {
    for (uint16_t j = 0; j < n; ++j)
    {
      printf("%9.3f ", mat[i * ldm + j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_im(const int* mat, uint16_t m, uint16_t n, uint16_t ldm, char* name)
{
  printf("%s = \n", name);
  for (uint16_t i = 0; i < m; ++i)
  {
    for (uint16_t j = 0; j < n; ++j)
    {
      printf("%d ", mat[i * ldm + j]);
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
      printf("%9.3f ", mat[j * ldm + i]);
    }
    printf("\n");
  }
  printf("\n");
}
