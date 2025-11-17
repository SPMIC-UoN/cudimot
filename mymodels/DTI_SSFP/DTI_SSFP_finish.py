#!/usr/bin/env python3

# A python code to convert the cudiMOT output from DTI_SSFP model into fsl-like dtifit outputParamaters to be passed from the DTI_SSFP output folder
#  In DTI model the parameters P are:
#  lam1: Param0
#  lam2: Param1
#  lam3: Param2
#  theta: Param3
#  phi:  Param4
#  psi:  Param5
#  S0:   Param6
#  Author: MKS mohamed.selim@nottingham.ac.uk"








import numpy as np
import nibabel as nib
import argparse
import os
import warnings


def load_nifti(file_path, mean_last_dim=False):
    img = nib.load(file_path)
    data = img.get_fdata()
    if mean_last_dim and data.ndim > 3:
        data = np.mean(data, axis=-1)
    return data, img.affine, img.header

def save_nifti(data, affine, header, out_path, intent_code=None):
    header = header.copy()
    threshold = 1e-6  # use 1e-6 or 1e-8 depending on your precision needs
    affine[np.abs(affine) < threshold] = 0
    if intent_code == 'vector':
        header.set_intent('vector', (), '')
    img = nib.Nifti1Image(data.astype(np.float32), affine, header)
    nib.save(img, out_path)


def compute_v2_components(S_theta, psi_data, v1x, v1y, v1z):
    warnings.filterwarnings("ignore")
    sqrt_S = np.sqrt(S_theta * S_theta)
    v2x_sn0 = (-np.cos(psi_data) * v1y - np.sin(psi_data) * v1x * v1z) / sqrt_S
    v2y_sn0 = (np.cos(psi_data) * v1x - np.sin(psi_data) * v1y * v1z) / sqrt_S
    v2z_sn0 = np.sin(psi_data) * sqrt_S

    v2x_s0 = -np.cos(psi_data)
    v2y_s0 = 0.0
    v2z_s0 = np.sin(psi_data)

    mask = (S_theta != 0).astype(int)
    inv_mask = (S_theta == 0).astype(int)

    v2x = v2x_sn0 * mask + v2x_s0 * inv_mask
    v2y = v2y_sn0 * mask + v2y_s0 * inv_mask
    v2z = v2z_sn0 * mask + v2z_s0 * inv_mask

    return v2x, v2y, v2z


class VectorGenerator:
    def __init__(self, vx, vy, vz, mask):
        self.data = np.zeros((*vx.shape[:3], 3))

        if vx.ndim == 3:
            self.data[..., 0] = vx * mask
            self.data[..., 1] = vy * mask
            self.data[..., 2] = vz * mask
        else:
            for k in range(vx.shape[0]):
                for l in range(vx.shape[1]):
                    for m in range(vx.shape[2]):
                        if mask[k, l, m] == 0.0:
                            continue
                        v_dyad = np.zeros((3, 3))
                        v_dyad[0, 0] = np.mean(vx[k, l, m, :] * vx[k, l, m, :])
                        v_dyad[0, 1] = np.mean(vx[k, l, m, :] * vy[k, l, m, :])
                        v_dyad[0, 2] = np.mean(vx[k, l, m, :] * vz[k, l, m, :])
                        v_dyad[1, 0] = np.mean(vy[k, l, m, :] * vx[k, l, m, :])
                        v_dyad[1, 1] = np.mean(vy[k, l, m, :] * vy[k, l, m, :])
                        v_dyad[1, 2] = np.mean(vy[k, l, m, :] * vz[k, l, m, :])
                        v_dyad[2, 0] = np.mean(vz[k, l, m, :] * vx[k, l, m, :])
                        v_dyad[2, 1] = np.mean(vz[k, l, m, :] * vy[k, l, m, :])
                        v_dyad[2, 2] = np.mean(vz[k, l, m, :] * vz[k, l, m, :])

                        v_dyad = np.nan_to_num(v_dyad, nan=0.0, posinf=0.0, neginf=0.0) ## added MKS
                        eigval, eigvec = np.linalg.eig(v_dyad)
                        self.data[k, l, m, :] = eigvec[:, np.argmax(eigval)]
                        #self.data = self.data * mask[..., None]  # Apply mask

def construct_tensor_from_spherical(lam1, lam2, lam3, v1, v2, v3):
    shape = lam1.shape
    D = np.zeros(shape + (3, 3))
    for i in range(3):
        for j in range(3):
            D[..., i, j] = (
                lam1 * v1.data[..., i] * v1.data[..., j] +
                lam2 * v2.data[..., i] * v2.data[..., j] +
                lam3 * v3[..., i] * v3[..., j]
            )
    return D


def compute_metrics(D, mask):
    shape = D.shape[:-2]
    D_flat = D.reshape(-1, 3, 3)
    evals = np.linalg.eigvalsh(D_flat)
    evals = np.clip(evals, 0, None)
    evals = evals.reshape(shape + (3,))

    lam1 = evals[..., 2]
    lam2 = evals[..., 1]
    lam3 = evals[..., 0]

    MD = (lam1 + lam2 + lam3) / 3
    FA = np.sqrt(1.5) * np.sqrt((lam1 - MD)**2 + (lam2 - MD)**2 + (lam3 - MD)**2) / \
         np.sqrt(lam1**2 + lam2**2 + lam3**2 + 1e-15)  # Added small constant to avoid division by zero
    RD = (lam2 + lam3) / 2
    AD = lam1

    return MD * mask, FA * mask, RD * mask, AD * mask, lam1 * mask, lam2 * mask, lam3 * mask


def main():
    parser = argparse.ArgumentParser(description="Compute DTI metrics for DWSSFP_DTI (FA, MD, AD, RD, V1, V2, V3, L1, L2, L3, S0) with mask.")
    parser.add_argument('--lam1', required=True)
    parser.add_argument('--lam2', required=True)
    parser.add_argument('--lam3', required=True)
    parser.add_argument('--theta', required=True)
    parser.add_argument('--phi', required=True)
    parser.add_argument('--psi', required=True)
    parser.add_argument('--s0', required=True)
    parser.add_argument('--mask', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    lam1, _, _ = load_nifti(args.lam1, mean_last_dim=True) # Load and mean over last dimension (50 volumes MCMC) if needed for lam1, lam2, lam3 and s0
    lam2, _, _ = load_nifti(args.lam2, mean_last_dim=True)
    lam3, _, _ = load_nifti(args.lam3, mean_last_dim=True)
    s0, _, _ = load_nifti(args.s0, mean_last_dim=True)
    theta, _, _ = load_nifti(args.theta)
    phi, _, _ = load_nifti(args.phi)
    psi, _, _ = load_nifti(args.psi)
    mask, aff, hdr = load_nifti(args.mask)
 
    # Compute V1
    S_theta = np.sin(theta)
    v1x = S_theta * np.cos(phi)
    v1y = S_theta * np.sin(phi)
    v1z = np.cos(theta)
    V1 = VectorGenerator(v1x, v1y, v1z, mask)

    # Compute V2
    v2x, v2y, v2z = compute_v2_components(S_theta, psi, v1x, v1y, v1z)
    V2 = VectorGenerator(v2x, v2y, v2z, mask)

    # Compute V3
    v3x = V1.data[..., 1] * V2.data[..., 2] - V2.data[..., 1] * V1.data[..., 2]
    v3y = V2.data[..., 0] * V1.data[..., 2] - V1.data[..., 0] * V2.data[..., 2]
    v3z = V1.data[..., 0] * V2.data[..., 1] - V2.data[..., 0] * V1.data[..., 1]
    V3 = np.stack([v3x, v3y, v3z], axis=-1) * mask[..., None]

    # Tensor and metrics
    D = construct_tensor_from_spherical(lam1, lam2, lam3, V1, V2, V3)
    MD, FA, RD, AD, L1, L2, L3 = compute_metrics(D, mask)

    # Save outputs
    save_nifti(FA, aff, hdr, os.path.join(args.outdir, 'dti_FA.nii.gz'))
    save_nifti(MD, aff, hdr, os.path.join(args.outdir, 'dti_MD.nii.gz'))
    save_nifti(AD, aff, hdr, os.path.join(args.outdir, 'dti_AD.nii.gz'))
    save_nifti(RD, aff, hdr, os.path.join(args.outdir, 'dti_RD.nii.gz'))
    save_nifti(L1, aff, hdr, os.path.join(args.outdir, 'dti_L1.nii.gz'))
    save_nifti(L2, aff, hdr, os.path.join(args.outdir, 'dti_L2.nii.gz'))
    save_nifti(L3, aff, hdr, os.path.join(args.outdir, 'dti_L3.nii.gz'))
    save_nifti(s0, aff, hdr, os.path.join(args.outdir, 'dti_S0.nii.gz')) # save the averaged S0 

    save_nifti(V1.data, aff, hdr, os.path.join(args.outdir, 'dti_V1.nii.gz'),)# intent_code='vector')
    save_nifti(V2.data, aff, hdr, os.path.join(args.outdir, 'dti_V2.nii.gz'),)# intent_code='vector')
    save_nifti(V3, aff, hdr, os.path.join(args.outdir, 'dti_V3.nii.gz'), )#intent_code='vector')


if __name__ == '__main__':
    main()
