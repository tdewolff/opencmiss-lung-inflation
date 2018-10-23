clc 
clear 
close all
addpath(genpath('/hpc/tdew803/Downloads/NIfTYPackages'));

img = 0 * ones(101, 101, 101, 'double');
nii = make_nii(img);
nii.hdr.dime.pixdim = [0 0.001 0.001 0.001 1 0 0 0];
save_nii(nii, 'cube_dy.nii');
nii = make_nii(img);
nii.hdr.dime.pixdim = [0 0.001 0.001 0.001 1 0 0 0];
save_nii(nii, 'cube_dz.nii');

for k = 1:101
    for j = 1:101
        for i = 1:101
            img(i,j,k) = (i-1)/100.0 * 0.02;
            if k==1 && j==1 && (i-1 == 0 || i-1 == 33 || i-1 == 66 || i-1 == 100)
                (i-1)/100 * 0.02
            end
        end
    end
end
nii = make_nii(img);
nii.hdr.dime.pixdim = [0 0.001 0.001 0.001 1 0 0 0];
save_nii(nii, 'cube_dx.nii');

figure;
imshow(squeeze(img(1,:,:)), []);