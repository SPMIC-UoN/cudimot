#define MACRO template <typename T> __device__ inline

MACRO void transpose_square_matrix(int N, T* A, T* C){
	for(int row=0;row<N;row++){
		for(int col=0;col<N;col++){
			C[row*N+col]=A[col*N+row];
		}
	}
}

MACRO void multiply_square_matrices(int N, T* A, T *B, T* C){
	for(int row=0;row<N;row++){
		for(int col=0;col<N;col++){
			C[row*N+col]=(T)0.0;
		}
	}
	for(int rowA=0;rowA<N;rowA++){
		for(int colB=0;colB<N;colB++){
			for(int elem=0;elem<N;elem++){
				C[rowA*N+colB]+= A[rowA*N+elem] * B[elem*N+colB];
			}
		}
	}
}
