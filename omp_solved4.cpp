/******************************************************************************
* FILE: omp_bug4.c
* DESCRIPTION:
*   This very simple program causes a segmentation fault.
* AUTHOR: Blaise Barney  01/09/04
* LAST REVISED: 04/06/05
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h> 
#define N 20				// COMMENT; The size of array that every thread had to handle was too large!! So makde N as 20. 
					//But in general, for N 1048 do the following :
					//	 (1) ulmit -s 100240 This is sets the stack size to a large value to stop overflow errors
					//       (2) export OMP_STACKSIZE=20m  This increases the OMP stack size 
					// 	*****BOTH these steps can be acheived by first sourcing the .sh file and then compiling and running this program*****
					// These two steps solves the SEG FAULT


int main (int argc, char *argv[]) 
{
	

int nthreads, tid, i, j;
double a[N][N];

/* Fork a team of threads with explicit variable scoping */
#pragma omp parallel shared(nthreads) private(i,j,tid,a)
  {

  /* Obtain/print thread info */
  tid = omp_get_thread_num();
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n", tid);

  /* Each thread works on its own private copy of the array */
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i][j] = tid + i + j;

  /* For confirmation */
  printf("Thread %d done. Last element= %f\n",tid,a[N-1][N-1]);

  }  /* All threads join master thread and disband */

}
