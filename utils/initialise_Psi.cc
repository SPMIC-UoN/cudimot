/*  initialise_Psi.cc

    This file initialise the fanning angle (psi) for fitting NODDI-Bingham model.
    
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

#define MyType double

using namespace std;

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
