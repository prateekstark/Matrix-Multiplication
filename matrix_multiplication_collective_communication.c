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
	if(argc == 2){
		int i, j, np, mype, ierr;
		double time_taken;
		struct timeval start, end;
		gettimeofday(&start, NULL); 
		int N = atoi(argv[1]);
		float *a = (float *)malloc((N*32) * sizeof(float *));
		float *b = (float *)malloc((32*N) * sizeof(float *));
		float *c = (float *)malloc((N*N) * sizeof(float *));


		ierr = MPI_Init(&argc, &argv);
		ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mype);

		if(mype == 0){
			
			for(i=0;i<N;i++){
				for(j=0;j<32;j++){
					// a[i*32+j] = drand48();
					a[i*32+j] = 0;
				}
			}

			for(i=0;i<32;i++){
				for(j=0;j<N;j++){
					// b[i*N+j] = drand48();
					b[i*N+j] = 0;
				}
			}
		}

		float *localCalcA = (float *)malloc((N*32/np) * sizeof(float *));
		float *localCRow = (float *)malloc((N*N/np) * sizeof(float *));
		MPI_Bcast(b, 32*N, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Scatter(a, N*32/np, MPI_FLOAT, localCalcA, N*32/np, MPI_FLOAT, 0, MPI_COMM_WORLD);

	
		Matrix_Multiply(localCalcA, b, localCRow, N/np, 32, N);
		MPI_Gather(localCRow, N*N/np, MPI_FLOAT, c, N*N/np, MPI_FLOAT, 0, MPI_COMM_WORLD);

		if(mype == 0){
			/*
			 // Validation Block

			float *c_check = (float *)malloc((N*N) * sizeof(float *));
			Matrix_Multiply(a, b, c_check, N, 32, N);
			
			printMatrix(c_check, N, N);

			printf("Is the result correct? %d\n", IsEqual(c, c_check, N, N));
			
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