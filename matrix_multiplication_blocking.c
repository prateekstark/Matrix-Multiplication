#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include<sys/time.h>

void Multiply_Serial(float *A, float *B, float *C, int m, int n, int p){
	int i, j, k;
	for (i = 0; i < m; i++){
		for (j = 0; j < p; j++){
			C[i*p + j] = 0;
			for (k = 0; k < n; k++)
				C[i*p + j] += A[i*n + k] * B[k*p + j];
		}
	}
}	

void Matrix_Multiply(float *A, float *B, float *C, int m, int n, int p, int np, int my_rank){
	int i,j,k;
	int work = m/np;
	if (my_rank == np-1 && (my_rank+1)*work < m) {
		for (i=my_rank*work; i<m; i++) {
			for (j = 0; j < p; j++){
				C[i*p + j] = 0;
				for (k = 0; k < n; k++)
					C[i*p + j] += A[i*n + k] * B[k*p + j];
			}
		}
	}
	else {
		for (i=my_rank*work; i<(my_rank+1)*work; i++){
			for (j = 0; j < p; j++){
				C[i*p + j] = 0;
				for (k = 0; k < n; k++)
					C[i*p + j] += A[i*n + k] * B[k*p + j];
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

int main(int argc, char *argv[]){
	if(argc == 2){
		int tag = 21;

		int N = atoi(argv[1]);
		float* A = (float*) malloc((N*32)*sizeof(float));
		float* B = (float*) malloc((32*N)*sizeof(float));
	    float* answer = (float*) malloc((N*N)*sizeof(float));
		float* C = (float*) malloc((N*N)*sizeof(float));

	    struct timeval start, end; 
		MPI_Init(NULL,NULL);
		int np;
		
		int my_rank;
		MPI_Comm_size(MPI_COMM_WORLD, &np);
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

		if(my_rank != 0) {

			MPI_Recv(&N,1,MPI_INT, 0, tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(A,N*32,MPI_FLOAT,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(B,32*N,MPI_FLOAT,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			int m=N, n = 32, p = N;
			
			for (int i=1;i<np;i++) {
				if (my_rank == i) {
					Matrix_Multiply(A,B,C,m,n,p,np,i);
					MPI_Send(C,(N*N),MPI_FLOAT,0,tag,MPI_COMM_WORLD);	
				}
			}
		}

		else{
			
			gettimeofday(&start, NULL); 

			for (int i=0; i<N*32; i++) {
				A[i] = (rand()/(float)RAND_MAX);
			}
			
			for (int i=0; i<32*N; i++) {
				B[i] = (rand()/(float)RAND_MAX);
			}
			
			for (int i=1;i<np;i++) {
				MPI_Send(&N,1,MPI_INT,i,tag,MPI_COMM_WORLD);
				MPI_Send(A,32*N,MPI_FLOAT,i,tag,MPI_COMM_WORLD);
				MPI_Send(B,32*N,MPI_FLOAT,i,tag,MPI_COMM_WORLD);
			}

			int m=N, n = 32, p = N;
			Matrix_Multiply(A,B,answer,m,n,p,np,my_rank);

			for (int i=1;i<np;i++){
				MPI_Recv(C,N*N,MPI_FLOAT,i,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				for (int j=0;j<N*N;j++) {
					if (C[j] != 0) {
						answer[j] = C[j];
					}
				}
			}

			gettimeofday(&end, NULL); 
	    	double time_taken;
	    	time_taken = (end.tv_sec - start.tv_sec) * 1e6; 
	    	time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
	        printf("Time Taken: %lf\n",time_taken);
	        float* answer_serial = (float*) malloc(N*N*sizeof(float));
			
			/* Validation Step
			// Multiply_Serial(A,B,answer_serial,N,32,N);
			
			// printMatrix(answer_serial, N, N);
			// printMatrix(answer, N, N);
			
			*/
			printf("Is my answer correct? %d\n", IsEqual(answer_serial, answer, N, N));
		}
		

		MPI_Finalize();

	}
	else{
		printf("Wrong Arguments!\n");
	}
	return 0;

}