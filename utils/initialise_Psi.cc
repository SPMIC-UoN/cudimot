/*  initialise_Psi.cc

    This file initialise the fanning angle (psi) for fitting NODDI-Bingham model.
    
    Moises Hernandez-Fernandez FMRIB Image Analysis Group
    
    Copyright (C) 2005 University of Oxford  */

#include <iostream>
#include <stdlib.h>
#include "newimage/newimageall.h"

#define MyType double

using namespace std;
using namespace NEWMAT;

int main (int argc, char* argv[]){
  cout << "Initialise Fanning angle - Psi" << endl;

  if(argc<5){
    cerr << "Usage: initialise_Psi Fibre_Direction(3D)  Fanning_Direction(3D) mask Output_Volume" << endl;
    exit(-1);
  }

  NEWIMAGE::volume4D<MyType> FibreDir;
  read_volume4D(FibreDir,argv[1]);

  NEWIMAGE::volume4D<MyType> FanningDir;
  read_volume4D(FanningDir,argv[2]);

  NEWIMAGE::volume<MyType> mask;
  read_volume(mask,argv[3]);

  if(FibreDir.tsize()!=3 || FanningDir.tsize()!=3){
    cerr << "Error: The input files must be 3D volumes" << endl;
  }

  NEWIMAGE::volume<MyType> Psi(FibreDir.xsize(),FibreDir.ysize(),FibreDir.zsize());
  Psi=0.0;
  
  Matrix mat(2,2);
  ColumnVector B(2);
  ColumnVector tmp(2);
  MyType theta,phi;

  for(int z=0;z<FibreDir.zsize();z++){
    for(int y=0;y<FibreDir.ysize();y++){
      for(int x=0;x<FibreDir.xsize();x++){

	if(mask(x,y,z)){

	  theta = acos(FibreDir[2](x,y,z));
	  phi = atan2(FibreDir[1](x,y,z), FibreDir[0](x,y,z));
	  
	  if(fabs(FibreDir[2](x,y,z))>0.1){
	    mat(1,1)=-sin(phi);
	    mat(1,2)=-cos(theta)*cos(phi); 
	    mat(2,1)=cos(phi);
	    mat(2,2)=-cos(theta)*sin(phi);
	    	    
	    B(1)=FanningDir[0](x,y,z);
	    B(2)=FanningDir[1](x,y,z);
	    
	    tmp = mat.i()*B;
	    
	    Psi(x,y,z) = atan2(tmp(2),tmp(1));
	    
	  }else if(fabs(sin(phi))>0.1){
	    mat(1,1)=-sin(phi);
	    mat(1,2)=-cos(theta)*cos(phi); 
	    mat(2,1)=0;
	    mat(2,2)=sin(phi);
	    	    
	    B(1)=FanningDir[0](x,y,z);
	    B(2)=FanningDir[2](x,y,z);
	    
	    tmp = mat.i()*B;
	    Psi(x,y,z) = atan2(tmp(2),tmp(1));
	    
	  }else{
	    mat(1,1)=cos(phi);
	    mat(1,2)=-cos(theta)*sin(phi); 
	    mat(2,1)=0;
	    mat(2,2)=sin(phi);
	    	    
	    B(1)=FanningDir[1](x,y,z);
	    B(2)=FanningDir[2](x,y,z);
	    
	    tmp = mat.i()*B;
	    Psi(x,y,z) = atan2(tmp(2),tmp(1));
	  }	  
	}	
      }
    }
  }

  string output;
  output.append(argv[4]);
    
  save_volume(Psi,output);
  
  cout << "Done" << endl;
  
  return 0;
}
