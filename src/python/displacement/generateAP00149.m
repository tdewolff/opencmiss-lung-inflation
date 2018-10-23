clc 
clear 
close all
addpath(genpath('/hpc/tdew803/Downloads/NIfTYPackages'));

img = 1 * ones(512, 512, 610, 'double');
nii = make_nii(img);
nii.hdr.dime.pixdim = [0 0.51 0.51 0.50 1 0 0 0];
save_nii(nii, 'AP00149_dx.nii');

img = 2 * ones(512, 512, 610, 'double');
nii = make_nii(img);
nii.hdr.dime.pixdim = [0 0.51 0.51 0.50 1 0 0 0];
save_nii(nii, 'AP00149_dy.nii');

img = 3 * ones(512, 512, 610, 'double');
nii = make_nii(img);
nii.hdr.dime.pixdim = [0 0.51 0.51 0.50 1 0 0 0];
save_nii(nii, 'AP00149_dz.nii');

%figure;
%imshow(squeeze(img(1,:,:)), []);