#!/usr/bin/env fslpython

############ convert dyads from cartesian to spherical (theta, phi)

# %% prep libs
import sys,os,argparse
import numpy as np
import math
import nibabel as nib
from tqdm import tqdm

# %% functions
def errchk(errflag):
    if errflag:
        print("Exit without doing anything..")
        quit()

def cart2sph(x,y,z):
    r =  np.sqrt(x*x + y*y + z*z)
    if r==0:
        theta = math.acos(z / 1)    # To avoid NaN when r==0
    else:
        theta = math.acos(z/r) 
    phi = math.atan2(y,x)  
    return theta, phi

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#### arguments
parser = MyParser(prog='dyadcart')
required = parser.add_argument_group('Required arguments')
required.add_argument("-dyad", metavar='<NIFTI>', dest='dyad_path', required=True)
required.add_argument("-mask", metavar='<NIFTI>', dest='mask_path', required=True)

argsa = parser.parse_args()
dyad_path = argsa.dyad_path
mask_path = argsa.mask_path

#### checks
errflag = 0
if not os.path.isfile(dyad_path):
    print(f"dyad file {dyad_path} not found")
    errflag = 1

if not os.path.isfile(mask_path):
    print(f"mask file {mask_path} not found")
    errflag = 1

errchk(errflag)

#### load dyad
dyad = nib.load(dyad_path)
dyad_img = dyad.get_fdata(dtype=np.float32)

#### load mask dyad
mask = nib.load(mask_path)
mask = mask.get_fdata(dtype=np.float32)

# Get indices where mask is 1
indices = np.argwhere(mask == 1)

# Prepare an array to hold the converted spherical coordinates
i, j, k, _ = dyad_img.shape 
sph_comps = np.zeros((i, j, k, 2), dtype=np.float32)

# Iterate over the relevant indices and convert each vector
for index in tqdm(indices, desc="Converting dyads to spherical coordinates"):
    i, j, k = index
    x = dyad_img[i, j, k, 0]
    y = dyad_img[i, j, k, 1]
    z = dyad_img[i, j, k, 2]
    v = cart2sph(x, y, z)
    sph_comps[i, j, k, :] = v

# Create new NIfTI images for the spherical coordinates (theta and phi) and save
theta_img = nib.Nifti1Image(sph_comps[:, :, :, 0], dyad.affine, dyad.header)
output_path = dyad_path.replace("dyads", "th")
nib.save(theta_img, output_path)

phi_img = nib.Nifti1Image(sph_comps[:, :, :, 1], dyad.affine, dyad.header)
output_path = dyad_path.replace("dyads", "ph")
nib.save(phi_img, output_path)

print('Done!')
