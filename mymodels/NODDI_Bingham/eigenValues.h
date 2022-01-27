#define MACRO template <typename T> __device__ inline

#define MPI 3.14159265358979323846
#define SQR(x) ((x)*(x)) // x^2

// calculate 3 EigenValues from a 3x3 Matrix.
// Matrix must be symmetric
// from wikipedia
// https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
// Smith, Oliver K. (April 1961), "Eigenvalues of a symmetric 3 × 3 matrix.", Communications of the ACM, 4 (4): 168, doi:10.1145/355578.366316

MACRO void getEigenValues_symm_3x3(T* A, T* EigenValues){
  // Given a real symmetric 3x3 matrix A, compute the eigenvalues
  T p1 = SQR(A[1]) + SQR(A[2]) + SQR(A[5]);
  if (p1 == (T)0.0){
    // A is diagonal.
    EigenValues[0] = A[0];
    EigenValues[1] = A[4];
    EigenValues[2] = A[8];
  }else{
    T q = (A[0]+A[4]+A[8])/(T)3.0; // trace(A)/3
    T p2 = SQR((A[0]-q)) + SQR((A[4]-q)) + SQR((A[8]-q)) + (T)2.0*p1;
    T p = sqrt_gpu(p2/(T)6.0);
    T B[9];  
    // B = (1 / p) * (A - q * I)       ... I is the identity matrix
    T p3 = (T)1.0/p; 
    B[0]=(A[0]-q)*p3;  B[1] = A[1]*p3;   	B[2] = A[2]*p3;
    B[3]=A[3]*p3;      B[4] = (A[4]-q)*p3; B[5] = A[5]*p3;
    B[6]=A[6]*p3;      B[7] = A[7]*p3;  	B[8] = (A[8]-q)*p3;
    
    //r = det(B) / 2
    T r = B[0]*((B[4]*B[8]) - (B[7]*B[5])) -B[1]*(B[3]*B[8] - B[6]*B[5]) + B[2]*(B[3]*B[7] - B[6]*B[4]);
    r = r/(T)2.0;
    
    // In exact arithmetic for a symmetric matrix  -1 <= r <= 1
    // but computation error can leave it slightly outside this range.
    T phi;
    if (r <= (T)-1.0){
      phi = (T)MPI / (T)3.0;
    }else if(r >= (T)1.0){
      phi = (T)0.0;
		 }else{
      phi = acos_gpu(r)/(T)3.0;
    }
    
    // the eigenvalues satisfy eig3 <= eig2 <= eig1
    EigenValues[0] = q + (T)2.0 * p * cos_gpu(phi);
    EigenValues[2] = q + (T)2.0 * p * cos_gpu(phi+((T)2.0*(T)MPI/(T)3.0));
    EigenValues[1] = (T)3.0 * q - EigenValues[0] - EigenValues[2];     // since trace(A) = eig1 + eig2 + eig3
  }
}
