#define MACRO template <typename T> __device__ inline

#define MPI 3.14159265358979323846
#define MSQRT3 1.73205080756887729352744634151
#define MINV3 0.33333333333333333
#define MINV54 0.018518518518519

MACRO T croot(T x){
  //if (x>=(T)0.0){ 
  //  return(pow_gpu(x,(T)MINV3));
  //}else{
  //  return -(pow_gpu(-x,(T)MINV3));
  //}
  if (x>=(T)0.0){ 
    return exp_gpu(log_gpu(x)*(T)MINV3);
  }else{
    return -exp_gpu(log_gpu(-x)*(T)MINV3);
  }
}

MACRO T find_t(T x1, T x2, T x3){
  T l1=-x1; 
  T l2=-x2; 
  T l3=-x3;

  T a3=l1*l2+l2*l3+l1*l3;
  T a2=(T)1.5-l1-l2-l3;
  T a1=a3-l1-l2-l3;
  T a0=(T)0.5*(a3-(T)2.0*l1*l2*l3);

  T p=(a1-a2*a2*(T)MINV3)*(T)MINV3;
  T q=(-(T)9.0*a2*a1+(T)27.0*a0+(T)2.0*a2*a2*a2)*(T)MINV54;
  T D=q*q+p*p*p;
  T offset=a2*(T)MINV3;

  T z1 = (T)0.0;
  T z2 = (T)0.0;
  T z3 = (T)0.0;
  if (D>(T)0.0){
    T ee=sqrt_gpu(D);
    T tmp=-q+ee; 
    z1=croot(tmp);
    tmp=-q-ee; 
    z1=z1+croot(tmp);
    z1=z1-offset; 
    z2=z1; 
    z3=z1;
  }else if(D<0){
    T ee=sqrt_gpu(-D);
    T tmp2=-q; 
    T angle=(T)2.0*MINV3*atan_gpu(ee/(sqrt_gpu(tmp2*tmp2+ee*ee)+tmp2));
    T tmp=cos_gpu(angle);
    tmp2=sin_gpu(angle);
    ee=sqrt_gpu(-p);
    z1=(T)2.0*ee*tmp-offset; 
    z2=-ee*(tmp+MSQRT3*tmp2)-offset; 
    z3=-ee*(tmp-MSQRT3*tmp2)-offset;
  }else{
    T tmp=-q;
    tmp=croot(tmp);
    z1=(T)2.0*tmp-offset; 
    z2=z1; 
    z3=z1;
    if (p!=(T)0.0 || q!=(T)0.0){
      z2=-tmp-offset; 
      z3=z2;
    }
  }
  z1=min_gpu(z1,z2); 
  z1=min_gpu(z1,z3);
  return z1;
}

MACRO T hyp_Sapprox(T x1, T x2, T x3){
  // Saddlepoint approximation of hypergeometric function of a matrix argument
  // Vector x has the eigenvalues 
  
  //sort(input,'descend');
  T tmp;
  if(x1>x2){
    if(x3>x1){
      tmp=x1;
      x1=x3;
      x3=tmp;
    }
  }else{ 
    if(x2>x3){
      tmp=x1;
      x1=x2;
      x2=tmp;
    }else{
      tmp=x1;
      x1=x3;
      x3=tmp;
    }
  }
  if(x3>x2){
    tmp=x2;
    x2=x3;
    x3=tmp;
  }

  T c1;  
  if (x1==(T)0.0 && x2==(T)0.0 && x3==(T)0.0){
    return 1;
  }else{
    T t=find_t(x1,x2,x3);
        
    T R=(T)1.0; 
    T K2=(T)0.0; 
    T K3=(T)0.0; 
    T K4=(T)0.0;
   
    R=R/sqrt_gpu(-x1-t);
    K2=K2+(T)1.0/((T)2.0*(x1+t)*(x1+t));
    K3=K3-(T)1.0/((x1+t)*(x1+t)*(x1+t));
    K4=K4+(T)3.0/((x1+t)*(x1+t)*(x1+t)*(x1+t));
       
    R=R/sqrt_gpu(-x2-t);
    K2=K2+(T)1.0/((T)2.0*(x2+t)*(x2+t));
    K3=K3-(T)1.0/((x2+t)*(x2+t)*(x2+t));
    K4=K4+(T)3.0/((x2+t)*(x2+t)*(x2+t)*(x2+t));
    
    R=R/sqrt_gpu(-x3-t);
    K2=K2+(T)1.0/((T)2.0*(x3+t)*(x3+t));
    K3=K3-(T)1.0/((x3+t)*(x3+t)*(x3+t));
    K4=K4+(T)3.0/((x3+t)*(x3+t)*(x3+t)*(x3+t));
    
    T T2=K4/((T)8.0*K2*K2)-(T)5.0*K3*K3/((T)24.0*K2*K2*K2);  
    c1=sqrt_gpu((T)2.0/K2)*MPI*R*exp_gpu(-t);
    T c3=c1*exp_gpu(T2);
    c1=c3/(T)4.0/(T)MPI;
  }
  
  return c1;
}
