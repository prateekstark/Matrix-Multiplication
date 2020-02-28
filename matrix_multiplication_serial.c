#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>


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

int main(int argc, char *argv[]){
	if(argc == 2){
		int i, j;
		double time_taken;
		struct timeval start, end;
		gettimeofday(&start, NULL); 
		int N = atoi(argv[1]);
		float *a = (float *)malloc((N*32) * sizeof(float *));
		float *b = (float *)malloc((32*N) * sizeof(float *));
		float *c = (float *)malloc((N*N) * sizeof(float *));

		for(i=0;i<N;i++){
			for(j=0;j<32;j++){
				a[i*32+j] = drand48();
			}
		}

		for(i=0;i<32;i++){
			for(j=0;j<N;j++){
				b[i*N+j] = drand48();
			}
		}
	
		Matrix_Multiply(a, b, c, N, 32, N);
		// printMatrix(c, N, N);
    

        gettimeofday(&end, NULL); 
		time_taken = (end.tv_sec - start.tv_sec) * 1e6;
		time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6; 

		printf("Time taken for Matrix Multiplication is: %lf\n", time_taken); 

	}
	else{
		printf("Wrong Arguments error!\n");
	}
	return 0;
}