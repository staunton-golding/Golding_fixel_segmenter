clear

%load data
wmfod_file=read_mrtrix('/Volumes/NO NAME/DWI/wm_3000t.mif');
wmfod_matlab=wmfod_file.data;

%use fixel direction transfrom from a lobe-based method (peak here)
direction_transform=read_mrtrix('/Volumes/NO NAME/DWI/pipeline/peak_fixels/directions.mif');
direction_transform_save=direction_transform.transform;

%mask to only segment fixels segmeneted with lobe-based method (for
%comparison, also selects WM voxels easily)
premask=read_mrtrix('/Volumes/NO NAME/DWI/pipeline/fixel.mif/index.mif');
premask=premask.data;

%load voxels that are used in WM RF estimation
voxels_matlab=read_mrtrix('/Volumes/NO NAME/DWI/rf_voxels_stride_corrected.mif');
voxels=voxels_matlab.data(:,:,:,3);

%initialize mask
mask=zeros(size(premask,1:3));

%fill in mask
mask(premask(:,:,:,1)~=0)=1;

%initialize dodecehedron, and tesselate 3 times - total of 1922 approximately
% equidistant sampled directions
TR=DodecahedronMesh; 
for i=1:3
	TR=SubdivideSphericalMesh(TR,1);
end
points=TR.Points;

%empty_mat is 1922x45 evaluation matrix of spherical harmonic expansion
%with all coefficients = 1 (empty_mat x coefficient vector gives 1x1922
%vector of SH expansion of any set of coefficient) (degree = 8, only even
%orders)
empty_mat=sh_by_matrix(points);

%model_odf is average ODFs of voxels used in WM RF estimation with only m=0
%coefficients kept (to preserve cylindrical symmetry)

%cap_thres is elevation wherein model_odf loses monotonicity (elevation
%with max x)

%single_fiber_max is max z value of model_odf
[model_odf,cap_thres,single_fiber_max]=single_fiber_odf(wmfod_matlab,voxels,points,empty_mat)

%run this to get mask of voxels adjacent to white matter, and the simulated
%fixels in these voxels. X_edge etc. is used to create box around white
%matter and adjacent voxels.
[vox_mask,x_edge,y_edge,z_edge,adj_fixels] = boundary_adjacent_wrapper(mask);

%find min and max x for box
x_start=min(x_edge);
x_end=max(x_edge);
clear x_edge
%y
y_start=min(y_edge);
y_end=max(y_edge);
clear y_edge
%z
z_start=min(z_edge);
z_end=max(z_edge);
clear z_edge

%segment fixels from voxels

%dir is an nx1 cell array where each cell input corresponds to mx4 matrix
%of values describing index and fixel information for each segmented voxel

%mask_total is a matrix combining true mask and simulated mask (used in 
%geometric fixel validation)

%total_fix is the total number of fixels found before validation

[dir, mask_total, total_fix] = run_fixel_wrapper(adj_fixels,mask,vox_mask,wmfod_matlab,model_odf,points,x_start, x_end, y_start, y_end, z_start, z_end,empty_mat,cap_thres,single_fiber_max);

%inputs into geometric validation
seperation_between_fixels=35;
number_to_compare=2;

%mag_fixel holds the magnitude of each found fixel

%dir_fixel holds the direction information of each fixel

%ind_fixel holds the indexing information of each fixel (starting fixel,
%number of fixels) (.mif file format)
[mag_fixel, dir_fixel, ind_fixel] = run_geometric_wrapper(dir,seperation_between_fixels, number_to_compare, x_start, x_end, y_start, y_end, z_start, z_end, mask_total);

%CREATE HISTOGRAMS

%load comparison fixel segmentations
mr_fix_SIFT=read_mrtrix('/Volumes/NO NAME/DWI/pipeline/SIFT_ss3t.mif/index.mif');
mr_fix_PEAK=read_mrtrix('/Volumes/NO NAME/DWI/pipeline/peak_fixels/index.mif');
mr_fix_PEAK=mr_fix_PEAK.data;
mr_fix_SIFT=mr_fix_SIFT.data;

%Find total number of fixels found with each modality
SIFT_fixels_found=sum(sum(sum(sum(mr_fix_comp_SIFT,4))))
PEAK_finding_fixels_found=sum(sum(sum(sum(mr_fix_comp_PEAK,4))))
my_found_fixels=sum(sum(sum(sum(ind_fixel,4))))

%voxels with at least 1 found fixel
idx_mine=ind_fixel>0;
idx_SIFT=mr_fix_SIFT>0;
idx_PEAK=mr_fix_PEAK>0;

%total percentage of WM voxels with >1 fixel
my_multi_fixel_percent=sum(sum(sum(ind_fixel(:,:,:,1)>1)))/size((find(mask)),1)*100
peak_multi_fixel_percent=sum(sum(sum(mr_fix_PEAK(:,:,:,1)>1)))/size((find(mask)),1)*100
SIFT_multi_fixel_percent=sum(sum(sum(mr_fix_SIFT(:,:,:,1)>1)))/size((find(mask)),1)*100

figure

histogram(comp_fix(idx_mine))
xlabel('Fixels per Voxel')
ylabel('Voxels')
title('my fixels')

histogram(mr_fix_comp_SIFT(idx_SIFT))
xlabel('Fixels per Voxel')
ylabel('Voxels')
title('SIFT fixels')

histogram(mr_fix_comp_PEAK(idx_PEAK))
xlabel('Fixels per Voxel')
ylabel('Voxels')
title('peak fixels')

%needed for python zero indexing (for mrview in MRTrix3)
ind_fixel(:,:,:,2)=ind_fixel(:,:,:,2)-1;
index.data=ind_fixel;

%transforms taken from original FODF and from lobe-based segmentation
index.vox=[2.5,2.5,2.5,1];
index.nfixels=int2str(my_found_fixels);
index.transform=wmfod_file.transform;

direc.data=dir_fixel;
direc.vox=[2.5,2.5,2.5,1];
direc.transform=direction_transform_save;

afd.data=mag_fixel;
afd.vox=[2.5,2.5,2.5,1];
afd.transform=direction_transform_save;

%write .mif files (insert name you want)
write_mrtrix(afd,'/Volumes/NO NAME/DWI/pipeline/test_fix_35_moving_50.mif/magnitude.mif');
write_mrtrix(direc,'/Volumes/NO NAME/DWI/pipeline/test_fix_35_moving_50.mif/directions.mif');
write_mrtrix(index,'/Volumes/NO NAME/DWI/pipeline/test_fix_35_moving_50.mif/index.mif');
%my nomenclature was test_fix_'seperation bewteen fixels'_'percent difference
%allowed (<)