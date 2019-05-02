function save_nii(ref_nii, data, output_nii)
% saves nifti files
% ref_nii is any refernce nifti file
% data is the volume/ data to be saved
% output_nii is the output file name
ref = nifti(ref_nii);
dat = file_array;
dat.fname = output_nii;
dimen=size(data);
dat.dim = dimen;
dat.dtype = 16;
out = nifti;
out.dat = dat;
out.mat = ref.mat;
out.mat_intent = ref.mat_intent;
out.mat0 = ref.mat0;
out.mat0_intent = ref.mat0_intent;
out.timing=ref.timing;
out.descrip='FSL5.0';
create(out);
out.dat(:,:,:,:) = data;
end

