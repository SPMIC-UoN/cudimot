/*  fanning2psi.cc

    This file get the psi angle from a 3D volume with Cartesian coordinates (x,y,z) indicating the Fanning-dispersion Oientation. 
    
    Moises Hernandez-Fernandez FMRIB Image Analysis Group
    
    Copyright (C) 2005 University of Oxford  */

#include <iostream>
#include <stdlib.h>
#include "newimage/newimageall.h"

#define M_PI 3.14159265358979323846

using namespace std;

int main (int argc, char* argv[]){
  cout << "Extracting angle Psi" << endl;
  
  if(argc<3){
    cerr << "Usage: fanning2psi mean_fibreOrientation(dyads) mean_fanningOrientation Output_basename" << endl;
    exit(-1);
  }

  NEWIMAGE::volume4D<double> fibreOr;
  read_volume4D(fibreOr,argv[1]);
  if(fibreOr.tsize()!=3){
    cerr << "Error: The input file must be 3D volumes (cartesian coordinates x,y,z)" << endl;
  }
  NEWIMAGE::volume4D<double> fanOr;
  read_volume4D(fanOr,argv[2]);
  if(fanOr.tsize()!=3){
    cerr << "Error: The input file must be 3D volumes (cartesian coordinates x,y,z)" << endl;
  }
  
  
  NEWIMAGE::volume<double> Psi(fibreOr.xsize(),fibreOr.ysize(),fibreOr.zsize());

  // Because there are several samples, probably fanning and fibre orientations are no orthogonal any more

  // Here we calculate 2 perpendical vectors and orthogonal to the fibre Orientation (x,y coordinates of fanning orientation on the ortoghonal plane)
  // Vector1 = dispersion (angle 0)
  // Vector2 = dispersion (angle pi/2)
  
  // Then x = dot(dispersionVector,Vector1)
  // y=dot(dispersionVector,Vector2)
  // psi=atan2(y,x)
	
  float th,ph,psi;
  ColumnVector cart;
  cart.ReSize(3);
  float Vector1[3];
  float Vector2[3];
  float Vx,Vy;
	
  for(int z=0;z<fibreOr.zsize();z++){
    for(int y=0;y<fibreOr.ysize();y++){
      for(int x=0;x<fibreOr.xsize();x++){
       
	cart(1)=fibreOr[0](x,y,z);
	cart(2)=fibreOr[1](x,y,z);
	cart(3)=fibreOr[2](x,y,z);
	
	cart2sph(cart,th,ph);
	
	psi=0;
	Vector1[0] = -cos(psi)*sin(ph) - cos(th)*cos(ph)*sin(psi);
	Vector1[1] = cos(ph)*cos(psi) - cos(th)*sin(ph)*sin(psi);
	Vector1[2] = sin(th)*sin(psi);
  
	psi=M_PI/2;
	Vector2[0] = -cos(psi)*sin(ph) - cos(th)*cos(ph)*sin(psi);
	Vector2[1] = cos(ph)*cos(psi) - cos(th)*sin(ph)*sin(psi);
	Vector2[2] = sin(th)*sin(psi);
     
	Vx=fanOr[0](x,y,z)*Vector1[0]+fanOr[1](x,y,z)*Vector1[1]+fanOr[2](x,y,z)*Vector1[2];
	Vy=fanOr[0](x,y,z)*Vector2[0]+fanOr[1](x,y,z)*Vector2[1]+fanOr[2](x,y,z)*Vector2[2];
	psi=atan2(Vy,Vx);

	Psi(x,y,z)=psi;	
      }
    }
  }

  string output1;
  output1.append(argv[3]);
  output1.append("_psi");
  
  save_volume(Psi,output1);

  cout << "Done" << endl;
  
  return 1;
}




