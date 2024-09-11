/*  cart2spherical.cc

    This file converts a 3D volume from Cartesian coordinates (x,y,z) to Spherical coordinates (th,ph). It uses the method cart2sph implemented in miscmaths.cc
    
    Moises Hernandez-Fernandez FMRIB Image Analysis Group
    
    Copyright (C) 2005 University of Oxford  */

#include <iostream>
#include <stdlib.h>
#include "newimage/newimageall.h"

using namespace std;
using NEWMAT::ColumnVector;
using MISCMATHS::cart2sph;

int main (int argc, char* argv[]){
  cout << "Cartesian coordinates to Spherical coordinates" << endl;
  
  if(argc<3){
    cerr << "Usage: cart2spherical Input_3D_volume Output_basename" << endl;
    exit(-1);
  }

  NEWIMAGE::volume4D<double> Cartesian;
  read_volume4D(Cartesian,argv[1]);

  if(Cartesian.tsize()!=3){
    cerr << "Error: The input file must be a 3D volume" << endl;
  }

  NEWIMAGE::volume<double> SphericalTH(Cartesian.xsize(),Cartesian.ysize(),Cartesian.zsize());
  NEWIMAGE::volume<double> SphericalPH(Cartesian.xsize(),Cartesian.ysize(),Cartesian.zsize());

  ColumnVector cart;
  cart.ReSize(3);
  float th,ph;

  for(int z=0;z<Cartesian.zsize();z++){
    for(int y=0;y<Cartesian.ysize();y++){
      for(int x=0;x<Cartesian.xsize();x++){
     
	cart(1)=Cartesian[0](x,y,z);
	cart(2)=Cartesian[1](x,y,z);
	cart(3)=Cartesian[2](x,y,z);
     
	cart2sph(cart,th,ph);

	SphericalTH(x,y,z)=th;
	SphericalPH(x,y,z)=ph;
	
      }
    }
  }

  string output1;
  output1.append(argv[2]);
  output1.append("_th");
  string output2;
  output2.append(argv[2]);
  output2.append("_ph");
  
  save_volume(SphericalTH,output1);
  save_volume(SphericalPH,output2);

  cout << "Done" << endl;
  
  return 1;
}
