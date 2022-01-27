#define MACRO template <typename T> __device__ inline

#define MPI 3.14159265358979323846
#define MSQRT3 1.73205080756887729352744634151
#define MINV3 0.33333333333333333
#define MINV54 0.018518518518519

MACRO T croot_doublePrecision(T x){
  if (x>=0.0){ 
    return exp(log(x)*MINV3);
  }else{
    return -exp(log(-x)*MINV3);
  }
}

MACRO double find_t_doublePrecision(T x1, T x2, T x3){
  double l1=-x1; 
  double l2=-x2; 
  double l3=-x3;

  double a3=l1*l2+l2*l3+l1*l3;
  double a2=1.5-l1-l2-l3;
  double a1=a3-l1-l2-l3;
  double a0=0.5*(a3-2.0*l1*l2*l3);

  double p=(a1-a2*a2*MINV3)*MINV3;
  double q=(-9.0*a2*a1+27.0*a0+2.0*a2*a2*a2)*MINV54;
  double D=q*q+p*p*p;
  double offset=a2*MINV3;

  double z1 = 0.0;
  double z2 = 0.0;
  double z3 = 0.0;
  if (D>0.0){
    double ee=sqrt(D);
    double tmp=-q+ee; 
    z1=croot_doublePrecision(tmp);
    tmp=-q-ee; 
    z1=z1+croot_doublePrecision(tmp);
    z1=z1-offset; 
    z2=z1; 
    z3=z1;
  }else if(D<0){
    double ee=sqrt(-D);
    double tmp2=-q; 
    double angle=2.0*MINV3*atan(ee/(sqrt(tmp2*tmp2+ee*ee)+tmp2));
    double tmp=cos(angle);
    tmp2=sin(angle);
    ee=sqrt(-p);
    z1=2.0*ee*tmp-offset; 
    z2=-ee*(tmp+MSQRT3*tmp2)-offset; 
    z3=-ee*(tmp-MSQRT3*tmp2)-offset;
  }else{
    double tmp=-q;
    tmp=croot_doublePrecision(tmp);
    z1=2.0*tmp-offset; 
    z2=z1; 
    z3=z1;
    if (p!=0.0 || q!=0.0){
      z2=-tmp-offset; 
      z3=z2;
    }
  }
  z1=min(z1,z2); 
  z1=min(z1,z3);
  return z1;
}

MACRO T hyp_Sapprox_doublePrecision(T x1, T x2, T x3){
  // Saddlepoint approximation of hypergeometric function of a matrix argument
  // Vector x has the eigenvalues 
  
  //sort(input,'descend');
  double tmp;
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

  double c1;  
  if (x1==0.0 && x2==0.0 && x3==0.0){
    return 1;
  }else{
    double t=find_t_doublePrecision(x1,x2,x3);

    double R=1.0; 
    double K2=0.0; 
    double K3=0.0; 
    double K4=0.0;

    R=R/sqrt(-x1-t);
    K2=K2+1.0/(2.0*(x1+t)*(x1+t));
    K3=K3-1.0/((x1+t)*(x1+t)*(x1+t));
    K4=K4+3.0/((x1+t)*(x1+t)*(x1+t)*(x1+t));
       
    R=R/sqrt(-x2-t);
    K2=K2+1.0/(2.0*(x2+t)*(x2+t));
    K3=K3-1.0/((x2+t)*(x2+t)*(x2+t));
    K4=K4+3.0/((x2+t)*(x2+t)*(x2+t)*(x2+t));
    
    R=R/sqrt(-x3-t);
    K2=K2+1.0/(2.0*(x3+t)*(x3+t));
    K3=K3-1.0/((x3+t)*(x3+t)*(x3+t));
    K4=K4+3.0/((x3+t)*(x3+t)*(x3+t)*(x3+t));
    
    double T2=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);  
    c1=sqrt(2.0/K2)*MPI*R*exp(-t);
    double c3=c1*exp(T2);
    c1=c3/4.0/MPI;
  }
  
  return c1;
}
