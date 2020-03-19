#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <iostream>
#include <fstream>
#include "utils.h"
#include <omp.h>


int main()
{   printf("DIMENSION     TIME \n");
    for(int n=10;n<=1000;n=n+10){
	//Timer t;

	//INITIALISING ALL VARIABLES ***********************************************
	int i,j,k		= 0;
	int N 			= n; // Dimension of Matrix
	double *f  		= (double *)malloc(N*N*sizeof(double)); // f_i matrix
	double *u1   		= (double *)malloc(N*N*sizeof(double)); // u_i matrix nth time step
	double *u2   		= (double *)malloc(N*N*sizeof(double)); // u_i matrix (n+1)th time step
	double *LHS   		= (double *)malloc(N*N*sizeof(double)); // -1*(delta u_i) matrix
	double *R   		= (double *)malloc(N*N*sizeof(double)); // R_i matrix to store residues
	int REPS    		= 15000;				// Repetitions/Iterations
	double h 		= 1.00/(N+1);				// Discretisation length
	double a,b,e,d		= 0.0;					// Grid point values
	double res,res_init 	= 0.0;					// Residue/Initial Residue
 
	//INITIALISING MATRICES **********************************************
	

	for(i=0;i<(N*N);i++){
		f[i]    = 1.0;
		u1[i]   = 0.0;
		u2[i]   = 0.0;
 		R[i]    = 0.0;
		LHS[i]  = 0.0;
	}


	//THE 2D JACOBI ALGORITHM $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	//t.tic();   
	double t = omp_get_wtime();    			// start timing

	for(int c=1;c<=REPS;c++){ 	// limit of # of iterations

		//#pragma omp parallel num_threads(10)
		//#pragma omp for collapse(2)
		for (int j=0;j<N;j++){
			for(int i=0;i<N;i++){

				if((i-1)<0){a = 0.0;}
				else { a =  u1[(i-1) + j*N] ; }

				if((j-1)<0){b = 0.0;}
				else { b =  u1[i + (j-1)*N] ; }

				if((i+1)==N){e = 0.0;}
				else { e =  u1[(i+1) + j*N] ; }

				if((j+1)==N){d = 0.0;}
				else { d =  u1[i + (j+1)*N] ; }

				u2[i + j*N] = 0.25 * (h*h*(f[i + j*N]) + a + b + e + d);
				a,b,e,d=0.0;
		
			}
			
		}
		
		// UPDATING ALL VALUES OF (n+1)th TIME STEP ################################

		#pragma omp parallel for schedule(static) num_threads(8)
		for(int i=0;i<(N*N);i++){
			u1[i] = u2[i];
			u2[i] = 0.0;
		}


		
		//COMPUTING RESIDUE ################################################
		
		#pragma omp parallel num_threads(8)
		#pragma omp for collapse(2)
		for (int j=0;j<N;j++){
		
			for(int i=0;i<N;i++){ 

				if((i-1)<0){a = 0.0;}
				else { a =  u1[(i-1) + j*N] ; }

				if((j-1)<0){b = 0.0;}
				else { b =  u1[i + (j-1)*N] ; }

				if((i+1)==N){e = 0.0;}
				else { e =  u1[(i+1) + j*N] ; }

				if((j+1)==N){d = 0.0;}
				else { d =  u1[i + (j+1)*N] ; }

				LHS[i + j*N] = (-1.00*(a + b + e + d) + 4.00*(u1[i + j*N]) )/(h*h);
			}
		
		}

		#pragma omp parallel for num_threads(8) reduction(+:res)		
		for (int i=0;i<(N*N);i++){
			R[i] = LHS[i] - f[i];
	       		res = res + (R[i])*(R[i]);
		
		} 
		res = pow((res),0.5);



		// CHECKING FOR CONVERGENCE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		if(c==1){res_init = res;}

		if(c>=1){
		
			if(floor(log10(res_init/res))==6){break;
			}
		}
		res = 0.0;	
	}

		// END OF JACOBI ITERATIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	printf("%d         %f \n",N,omp_get_wtime()-t); //t.toc()); 

	// FREE ALL ALLOCATIONS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	free(u1);
	free(u2);
	free(f);
	free(LHS);
	free(R);
   }
	return 0;
}
