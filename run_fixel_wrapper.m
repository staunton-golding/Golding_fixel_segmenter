function [dir,mask_total, total_fix, number_of_fifty_uses] = run_fixel_wrapper(adj_fixel,mask,vox_mask,wmfod_matlab,model_odf,points, x_start, x_end, y_start, y_end, z_start, z_end,empty_mat,cap_thres,single_fiber_max)
% Usage: [dir,mask_total] = run_fixel_wrapper(adj_fixel,mask,vox_mask,wmfod_matlab,model_odf,points, x_start, x_end, y_start, y_end, z_start, z_end,empty_mat,cap_thres,single_fiber_max)
%
% run_fixel_wrapper iterates the function 'fixel_wrapper' over the
% entirety of the whole brain white matter (or simulated) mask, and outputs
% from this function are stored in a data structure easily able to be converted to
% the .mif fixel data format (that conversion is done in a later function (run_geometric_wrapper)).
% Many of the inputs are the same as the function 'fixel_wrapper'
%
% Input Variables
%
%   
%   adj_fixels: nx3 matrix of xyz directions of each simulated fixel where
%               the nth row corresponds to the nth written simulated fixel
%
%   mask: binary mask of the DMRI scan, where values = +-1 correspond to the
%         white matter / area to be segmented, and all values = 0
%         correspond to areas to be ignored
%
%   vox_mask: mask of DMRI scan where values = -2 correspond voxels with
%             simulated fixels. -2 was chosen to differentiate simulated
%             fixels / voxels from 'real' fixels / voxels.
%
%   wmfod_matlab: 4D matrix where first three dimensions correspond to xyz
%                 indices and 4th dimension corresponds to SH coeffecients
%                 of FODFs (45 coefficients per FODF)
%   
%   model_odf: SH coefficient representation of single-fiber model ODF (all
%              non m=0 coefficients should = 0, to maintain cylindrical symmetry
%
%   points: an nx3 matrix of sampled points around a unit sphere
%
%   x_start: smallest x_value in vox_mask (all real voxels have x idx >=)
% 
%   x_end: largest x_value in vox_mask (all real voxels have x idx <=)
%
%   y_start: smallest y_value in vox_mask (all real voxels have y idx >=)
%
%   y_end: largest y_value in vox_mask (all real voxels have y idx <=)
%
%   z_start: smallest z_value in vox_mask (all real voxels have x idx >=)
%   
%   z_end: largest z_value in vox_mask (all real voxels have x idx <=)
%
%   empty_mat: nx45 matrix where index (n,m) correlates to the nth spherical
%              harmonic expansion at point m, where each column is an
%              evaluated spherical harmonic, stored in 'standard' order:
%              (l,m) = (0,0), (1,-1), (1,0), (1,1), (2,-1), ...
%
%   cap_thres: The elevation wherein, if starting at max radial value, the
%              single-fiber model FODF loses its monotonicity
%
%   single_fiber_max: the max radial value of the single-fiber model ODF
%
% Output Variables
%
%   dir: an nx1 cell array where each cell input corresponds to mx4 matrix
%        of values describing index and fixel information for each segmented
%        voxel
%
%   mask_total: a matrix combining true mask and simulated mask (used in
%               geometric fixel validation)
%
%   total_fix: the total number of fixels found before validation
%
% note* only x_start, x_end, y_start, y_end, z_start, z_end, points,
% wmfod_matlab, and mask AREN'T direct outputs of previously called functions



%combine masks together
mask_total=mask+vox_mask;

%initialize cell array for all fixels (simulated and not simulated)
dir={};

%counter
dir_idx=0;

%the adjacent fixels should be encountered in the order in which they were
%originally written given consistent strides, need counter to access
%correct index in array of simulated fixels.
adj_idx=0;
total_fix=0;
fifty_percent_uses=0;
%start and stop loop indices form a box around the area of the scan known
%to contain fixels and simulated fixels
for jj=x_start:x_end
    for kk=y_start:y_end
        for ll=z_start:z_end
            
            
            %increase counter
            dir_idx=dir_idx+1;

            %record index of current voxel of interest
            idx=[jj kk ll];
            
            %if voxel is in the white matter, call fixel_wrapper function
            if abs(mask(jj,kk,ll))==1
                test_odf=wmfod_matlab(jj,kk,ll,:);
                test_odf=squeeze(test_odf);
                [directions_scaling, ~,fifty_one_voxel]=fixel_wrapper(points,test_odf,model_odf,empty_mat,cap_thres,single_fiber_max);
                fifty_percent_uses= fifty_percent_uses+fifty_one_voxel;
                number_of_fixels=size(directions_scaling,1);

            %if voxel contains simulated fixels (boundary voxels), record
            %simulated fixel info into cell array
            elseif mask_total(jj,kk,ll)==-2
                adj_idx=adj_idx+1;
                directions_scaling=[adj_fixel(adj_idx,:),1];
                %so that simulated voxels can be identified
                number_of_fixels=-2;
            
            %in non white matter and non adjacent voxels, just record
            %number of fixels and zero and directions/scaling as empty
            else
                directions_scaling=[];
                number_of_fixels=0;
            end
            
            %Each cell input is recorded as first row = x y z and number of fixels then
            %lower rows as fixel direction + scaling factor
            if number_of_fixels>0
            total_fix=total_fix+number_of_fixels;
            end
            dir{dir_idx}=[idx,number_of_fixels;directions_scaling];
        end
    end
end

%how many fixels are found with poorly fitting caps while still being above
%size threshold
number_of_fifty_uses=fifty_percent_uses;

end
