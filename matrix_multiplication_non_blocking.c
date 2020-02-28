#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<mpi.h>

// We have used pointers to denote matrix as, when we use matrix, for large N, it gave segmentaton fault.

void printMatrix(float *A, int m, int n){
	int i, j;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			printf("%f ", A[i*n + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void Matrix_Multiply(float *A, float *B, float *C, int m, int n, int p){
	int i, j, k;
	for(i=0; i<m; i++){
		for(j=0; j<p; j++){
			C[i*p + j] = 0;
			for(k=0; k<n; k++){
				C[i*p + j] += A[i*n + k]*B[k*p + j];
			}
		}
	}
}

int IsEqual(float *A, float *B, int m, int n){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			if(A[i*n + j] != B[i*n + j]){
				printf("%f\n", A[i*n + j]);
				printf("%f\n", B[i*n + j]);
				return 0;
			}
		}
	}
	return 1;
}

void Transpose(float *A, float* B, int m, int n){
	int i, j;
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			B[j*m+i] = A[i*n+j];
		}
	}
}

int main(int argc, char *argv[]){
	// printf("Hey there! I just got started.\n");

	if(argc == 2){
		int i, j, np, mype, ierr;
		double time_taken;
		struct timeval start, end;
		gettimeofday(&start, NULL); 
		int N = atoi(argv[1]);
		float *a = (float *)malloc((N*32) * sizeof(float *));
		float *b = (float *)malloc((32*N) * sizeof(float *));
		float *b_column = (float *)malloc((32*N) * sizeof(float *));
		float *c = (float *)malloc((N*N) * sizeof(float *));

		ierr = MPI_Init(&argc, &argv);
		ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mype);
		
		MPI_Request request;
		MPI_Status status;
		// double start_time = MPI_Wtime();
		
		if(mype == 0){
			
			for(i=0;i<N;i++){
				for(j=0;j<32;j++){
					a[i*32+j] = drand48();
				}
			}

			for(i=0;i<32;i++){
				for(j=0;j<N;j++){
					b[i*N+j] = drand48();
					b_column[j*32+i] = b[i*N+j];
				}
			}
		}
		
		//Now we will begin scattering, dividing the n/p blocks among all processes

		float localCalcA[N*32/np];
		float localCalcB[32*N/np];

		MPI_Scatter(a, N*32/np, MPI_FLOAT, localCalcA, N*32/np, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Scatter(b_column, 32*N/np, MPI_FLOAT, localCalcB, 32*N/np, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		j = mype;
		int k = 0;
		float sum = 0;

		float localCalcC[N*N/np];

		for(int iter = 0; iter<np; iter++){
			for(i=0; i<N/np;i++){
				for(k=0;k<N/np;k++){
					sum = 0;
					for(int f=0;f<32;f++){
						sum += localCalcA[f + i*32]*localCalcB[f+k*32];
					}
					localCalcC[k*N+i+j*N/np] = sum; 
				}
			}
			if(j == np-1){
				j = 0;
			}
			else{
				j++;
			}
			MPI_Isend(localCalcA, N*32/np, MPI_FLOAT, (mype==0 ? np-1:mype-1), 0, MPI_COMM_WORLD, &request);
			MPI_Irecv(localCalcA, N*32/np, MPI_FLOAT, (mype==np-1 ? 0:mype+1), 0, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
		}

		MPI_Gather(localCalcC, N*N/np, MPI_FLOAT, c, N*N/np, MPI_FLOAT, 0, MPI_COMM_WORLD);
		if(mype == 0){
			/* Validation Block

			float *c_check = (float *)malloc((N*N) * sizeof(float *));
			float *c_transpose = (float *)malloc((N*N) * sizeof(float *));
			Transpose(c, c_transpose, N, N);
			
			Matrix_Multiply(a, b, c_check, N, 32, N);
			
			printMatrix(c_transpose, N, N);
			printMatrix(c_check, N, N);

			printf("Is the result correct? %d\n", IsEqual(c_transpose, c_check, N, N));
			
			*/

	    	gettimeofday(&end, NULL); 
			time_taken = (end.tv_sec - start.tv_sec) * 1e6;
			time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6; 
			printf("Time taken for Matrix Multiplication is: %lf\n", time_taken); 
		}

    	MPI_Finalize();

	}
	else{
		printf("Wrong Arguments error!\n");
	}
	return 0;
}