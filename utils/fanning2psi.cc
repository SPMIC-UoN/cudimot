/*  fanning2psi.cc

    This file get the psi angle from a 3D volume with Cartesian coordinates (x,y,z) indicating the Fanning-dispersion Oientation. 
    
    Moises Hernandez-Fernandez FMRIB Image Analysis Group
    
    Copyright (C) 2005 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk

    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK


    LICENCE

    FMRIB Software Library, Release 6.0 (c) 2018, The University of
    Oxford (the "Software")

    The Software remains the property of the Oxford University Innovation
    ("the University").

    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.

    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.

    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.

    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    fsl@innovation.ox.ac.uk quoting Reference Project 9564, FSL.*/

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




