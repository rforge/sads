#include <R.h>

void quadrado(int *n)
{
  int i, k = *n;

  k = k*k;
  
  Rprintf("%d\n", k);
}

