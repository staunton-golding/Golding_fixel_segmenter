function [vox_mask,x_edge,y_edge,z_edge,adj_fixels]=boundary_adjacent_wrapper(mask)
% Usage: [vox_mask,x_edge,y_edge,z_edge,adj_fixels]=boundary_adjacent_wrapper(mask)
%
% boundary_adjacent_wrapper places simulated fixels in voxels adjacent to
% the white matter (or simulated) mask. These simulated voxels are placed
% normal to the white matter, and are used for validation of fixels found
% in voxels on the boundary of the white matter / mask. Fixels are believed
% to enter the gray matter in a direction perpendicular to the gray matter.
%
% Input Variable:
%   
%   mask: binary mask of the DMRI scan, where values = +-1 correspond to the
%         white matter / area to be segmented, and all values = 0
%         correspond to areas to be ignored
%
% Output Variables:
%
%   vox_mask: mask of DMRI scan where values = -2 correspond voxels with
%             simulated fixels. -2 was chosen to differentiate simulated
%             fixels / voxels from 'real' fixels / voxels.
%   x_edge: x_values corresponding to voxels in vox_mask
%   y_edge: y_values corresponding to voxels in vox_mask
%   z_edge: z_values corresponding to voxels in vox_mask
%   adj_fixels: nx3 matrix of xyz directions of each simulated fixel where
%               the nth row corresponds to the nth written simulated fixel



%take gradient of mask (populated in voxels adjacent to white matter) 
%this produces fixels normal to surface in voxels adjacent to white matter
[gmag,gazimuth,gelevation]=imgradient3(mask,'prewitt');

%make gmag a consistent magnitude = 1 so that the fixels are represented as
%unit vectors
gmag((gmag~=0))=1;
gmag=gmag-mask;
%use gmag also to find boundary voxels with fixels (to ignore interior voxels)
gmag(gmag<0)=0;


%reshape all three into 1d arrays
gazimuth=reshape(gazimuth,[1,96*96*60]);
gelevation=reshape(gelevation,[1,96*96*60]);
gmag=reshape(gmag,[1,96*96*60]);

%identify boundary voxels with simulated fixels
vox_fix_idx=find(gmag);
vox_mask=zeros(size(mask));

%so that the simulated voxels can be identified as simulated
vox_mask(vox_fix_idx)=-2;

%get x, y, z directions for each simulated fixel
[x, y, z]=sph2cart(gazimuth(vox_fix_idx),gelevation(vox_fix_idx),gmag(vox_fix_idx));
adj_fixels=[x',y',z'];

%for use in constructing box around white matter + adjacent voxels
[x_edge,y_edge,z_edge]=ind2sub([96,96,60],vox_fix_idx);

end