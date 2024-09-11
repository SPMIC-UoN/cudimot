/*  Copyright (C) 2004 University of Oxford  */

#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;
using namespace NEWMAT;


int main ( int argc, char **argv ){
  if(argc<6){
    cout<<"usage: getFanningOrientation <theta_vol> <phi_vol> <psi_vol> <mask> <output>"<<endl;
    exit(1);
  }

  volume4D<float> ths,phs,psis;
  read_volume4D(ths,argv[1]);
  read_volume4D(phs,argv[2]);
  read_volume4D(psis,argv[3]);
  volume<float> mask;
  string oname;
  
  read_volume(mask,argv[4]);
  oname=argv[5];
  
  volume4D<float> fanning(ths.xsize(),ths.ysize(),ths.zsize(),3);
  fanning=0;
  copybasicproperties(ths,fanning);


  SymmetricMatrix dyad(3);dyad=0;
  ColumnVector dir(3);
  
  DiagonalMatrix dyad_D; //eigenvalues
  Matrix dyad_V; //eigenvectors

  for(int k=ths.minz();k<=ths.maxz();k++){
    for(int j=ths.miny();j<=ths.maxy();j++){
      for(int i=ths.minx();i<=ths.maxx();i++){
	if(mask(i,j,k)!=0){

	  dyad=0;
	  for(int s=ths.mint();s<=ths.maxt();s++){
	    float th=ths(i,j,k,s);
	    float ph=phs(i,j,k,s);
	    float psi=psis(i,j,k,s);
	    dir(1)=-cos(psi)*sin(ph) - cos(th)*cos(ph)*sin(psi); 
	    dir(2)=cos(ph)*cos(psi) - cos(th)*sin(ph)*sin(psi);
	    dir(3)=sin(th)*sin(psi);
	    dyad << dyad+dir*dir.t();
	  }
	  dyad = dyad/float(ths.maxt()-ths.mint()+1);
	  EigenValues(dyad,dyad_D,dyad_V);
	  int maxeig;
	  if(dyad_D(1)>dyad_D(2)){
	    if(dyad_D(1)>dyad_D(3)) maxeig=1;
	    else maxeig=3;
	  }
	  else{
	    if(dyad_D(2)>dyad_D(3)) maxeig=2;
	    else maxeig=3;
	  }
	  fanning(i,j,k,0)=dyad_V(1,maxeig);
	  fanning(i,j,k,1)=dyad_V(2,maxeig);
	  fanning(i,j,k,2)=dyad_V(3,maxeig);
	}else{
	  fanning(i,j,k,0)=0;
	  fanning(i,j,k,1)=0;
	  fanning(i,j,k,2)=0;
	}
      }
    }
  }

  fanning.setDisplayMaximumMinimum(1,-1);
  save_volume4D(fanning,oname);
  
  return 0;
}









