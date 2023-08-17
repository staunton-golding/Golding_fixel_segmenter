function [mag_fixel, dir_fixel, ind_fixel] = run_geometric_wrapper(dir,seperation_between_fixels, number_to_compare, x_start, x_end, y_start, y_end, z_start, z_end, mask_total)
% Usage: [mag_fixel, dir_fixel, ind_fixel] = run_geometric_wrapper(dir,seperation_between_fixels, number_to_compare, x_start, x_end, y_start, y_end, z_start, z_end, mask_total)
%
% run_geometric_wrapper sifts through all found segmented fixels, and
% validates fixels that are shown to be likely to exist via geometric
% validation. Our definition for geometric validation is that a fixel in a
% voxel should have a corresponding fixel in at least two adjacent voxels.
% The number of coherently oriented adjacent fixels needed for validation
% as well as the angle from which a fixel is considered to be coherently
% oriented is user-defined. The function outputs the data portion for the
% direction, magnitude, and index .mif files to be written
%
% To write a .mif file without any filtering, simply set number_to_compare = 0 
%
% Input Variables:
%
%   dir: an nx1 cell array where each cell input corresponds to mx4 matrix
%        of values describing index and fixel information for each segmented
%        voxels
%
%   seperation_between_fixels: The angle between the fixel being validated
%                              and adjacent fixels for the two fixels to be
%                              considered coherently oriented (how much can
%                              a fixel change direction in the span of one
%                              voxel)
%
%   number_to_compare: How many fixels must be coherently oriented to the
%                      fixel being validated for that target fixel to be
%                      considered likely
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
%   mask_total: a matrix combining true mask and simulated mask (used in
%               geometric fixel validation)
%
% Ouptut Variables:
%
%   mag_fixel: 1xn array where the nth input corresponds to the nth written
%              fixel's magnitude (derived from scaling needed for fit of single fiber ODF)
%
%   dir_fixel: 3xn array where nth column corresponds to nth fixel's
%              direction (first row = x, second = y, third = z).
%
%   ind_fixel: the index data for the sparse representation of fixel data.
%              A 4d matrix where the first three dimensions correspond to
%              xyz indices and the 4th dimension has two inputs. The first
%              input is the number of fixels in the voxel, the second input
%              is the starting fixel for that voxel (starting fixel in the
%              mag_fixel and dir_fixel arrays).

%validation is run through twice (in case a fixel is used for validation
%only to be deleted later
for kk=1:2

%loop through 'box'
for idx_x=x_start:x_end
    for idx_y=y_start:y_end
        for idx_z=z_start:z_end
            voxel_of_interest=[idx_x,idx_y,idx_z];
            
            %only validate non simulated fixels
            if abs(mask_total(idx_x,idx_y,idx_z))==1
                %geometric validation
                out_fixel=geometric_wrapper(dir,voxel_of_interest,seperation_between_fixels, number_to_compare, x_start, x_end, y_start, y_end, z_start, z_end);
          
                %update cell array each time fixel is checked (if a fixel is
                %deleted, update cell array with that information)
                
                if size(out_fixel,1)>1
                    dir{sub2ind([z_end-z_start+1,y_end-y_start+1,x_end-x_start+1],idx_z-z_start+1,idx_y-y_start+1,idx_x-x_start+1)}=out_fixel;
                else
                    dir{sub2ind([z_end-z_start+1,y_end-y_start+1,x_end-x_start+1],idx_z-z_start+1,idx_y-y_start+1,idx_x-x_start+1)}=[0,0,0,0];
                end
            
            end
           
        end
    end
end
end

%empty matrices for mif output files
dir_fixel=[];
mag_fixel=[];

%zero array for 4d .mif fixel output
ind_fixel=zeros(96,96,60,2);

%start counting of fixels 1 (for python indexing starts at 0 - this is
%corrected later in pipeline.m)
start_fix=1;

for mm=1:1:length(dir)

    %access the fixel info in each voxel
    temp_fix=dir{mm};

    %ignore simulated fixels
    if temp_fix(1,4)>0
  
        %in case a fixel has been deleted, count number of fixels in each
        %voxel this way (instead of using temp_fix(1,4)) 
        number_fix=size(temp_fix,1)-1;
        
        if mm==1
            %update the direction.mif file with each new fixel (sparse data rep)
            dir_fixel=temp_fix(2:end,1:3)';
        
            %update the magnitude.mif file with each new fixel (sparse data rep)
            mag_fixel=temp_fix(2:end,4)';
        else
            %update the direction.mif file with each new fixel (sparse data rep)
            dir_fixel=[dir_fixel,temp_fix(2:end,1:3)'];
            
            %update the magnitude.mif file with each new fixel (sparse data rep)
            mag_fixel=[mag_fixel,temp_fix(2:end,4)'];
        end
    
    %get the xyz coordinates of the voxel (inside the box)

    [z_fix,y_fix,x_fix]=ind2sub([z_end-z_start+1,y_end-y_start+1,x_end-x_start+1],mm);

    %record fixel info in 4th dimension
    %xyz has to be shifted in accordance with box
    ind_fixel(x_start+x_fix-1,y_start+y_fix-1,z_start+z_fix-1,:)=[number_fix,start_fix]; %-1 bc box start at fixel containing voxel??

    %update starting point for next voxel
    start_fix=start_fix+number_fix;
    end
 
    clear temp_fix
end

%transpose matrices (for .mif format cohesion)
dir_fixel=dir_fixel';
mag_fixel=mag_fixel';

end
