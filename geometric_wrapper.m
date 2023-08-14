function [output_fixel] = geometric_wrapper(dir,voxel_of_interest,seperation_between_voxels, number_to_compare, x_start, x_end, y_start, y_end, z_start, z_end)
% Usage: [output_fixel] = geometric_wrapper(dir,voxel_of_interest,seperation_between_voxels, number_to_compare, x_start, x_end, y_start, y_end, z_start, z_end)
%
% geometric_wrapper2 validates fixels that are shown to be likely to exist via geometric
% validation. Our definition for geometric validation is that a fixel in a
% voxel should have a corresponding fixel in at least two adjacent voxels.
% The number of coherently oriented adjacent fixels needed for validation
% as well as the angle from which a fixel is considered to be coherently
% oriented is user-defined. 
%
% Input Variables:
%
%   dir: an nx1 cell array where each cell input corresponds to mx4 matrix
%        of values describing index and fixel information for each segmented
%        voxels
%
%   voxel_of_interest: the xyz index of the voxel whose fixels are being
%                      validated
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
% Output Variable:
%
%   output_fixel: A nx4 matrix representing n-1 fixels. The first row's
%                 values correlate to [x_index, y_index, z_index, number of
%                 fixels] while each following row's values correlate to
%                 [x_dir, y_dir, z_dir, scaling factor] for each validated
%                 fixel.
    
%identify xyz idx of voxel of interest
x_vox=voxel_of_interest(1);
y_vox=voxel_of_interest(2);
z_vox=voxel_of_interest(3);

%plus 1 needed bc first index would be x_start-x_vox which = 0 (when
%x_start=v_vox) but needs to equal 1
main_fix_cart=squeeze(dir{sub2ind([z_end-z_start+1,y_end-y_start+1,x_end-x_start+1],z_vox-z_start+1,y_vox-y_start+1,x_vox-x_start+1)});
main_fix_cart=main_fix_cart(1:end,1:4);

focus_fix=0;

%loop through all fixels in voxel of interest
for mm=1:(size(main_fix_cart,1)-1)

focus_fix=focus_fix+1;

%matrix to record angles of seperation between fixels
fix_ang_sep_mat=zeros(27,27);

%counter
jj=0;
%loop through indices of surrounding voxels
for idx_1=x_vox-1:x_vox+1

    for idx_2=y_vox-1:y_vox+1
        
        for idx_3=z_vox-1:z_vox+1
            
            %get temporary voxel of fixels (fixels in neighboring voxel)
            fix_temp=dir{sub2ind([z_end-z_start+1,y_end-y_start+1,x_end-x_start+1],idx_3-z_start+1,idx_2-y_start+1,idx_1-x_start+1)};
            %increase counter
            jj=jj+1;

            %coordinate of the surrounding voxels (centered around main
            %voxel)
            vox_coord=[idx_1-x_vox,idx_2-y_vox,idx_3-z_vox];
         
            %angle from considered fixel to surrounding voxels
            vec1=main_fix_cart(mm+1,1:3)';
            vec2=vox_coord';
            ang_from_vox=atan2d(norm(cross(vec1,vec2)),dot(vec1,vec2));

            %number of fixels in neighboring voxel
            num_fix_temp=size(fix_temp(2:end,1:3),1);
            %array to hold angles from considered fixel to fixels of 1 adjacent voxel
            temp_angle_from_mid=zeros(1,num_fix_temp);
 
            %angles between fixels
            vec_main_fix=main_fix_cart(mm+1,1:3);
            for uu=1:num_fix_temp
               vec_temp_fix=fix_temp(uu+1,1:3);
               temp_angle_from_mid(uu) = atan2d(norm(cross(vec_main_fix,vec_temp_fix)),dot(vec_main_fix,vec_temp_fix));
            end

            %normalize to equal -90 to 90 deg (not 0 to 180)
            if ang_from_vox>90
                ang_from_vox=ang_from_vox-180;
            end
           
            %concatenate
            %normalize to equal -90 to 90 deg (not 0 to 180)
            temp_angle_from_mid(temp_angle_from_mid>90)=180-temp_angle_from_mid(temp_angle_from_mid>90);
            angles=[ang_from_vox, temp_angle_from_mid];
            %first 3 entries are coordinate of adjacent voxel
            fix_ang_sep_mat(jj,1:3)=[idx_1,idx_2,idx_3];
            %next are angles between fixel of interest and fixels in 1
            %adjacent voxel
            fix_ang_sep_mat(jj,4:length(angles)+3)=angles;
            fix_ang_sep_mat(jj,length(angles)+4:end)=NaN;

            clear fix_temp
            clear temp_angle_from_mid
          
        end
    end
end


%ignore same-voxel info for validation
fix_ang_sep_mat(14,:)=[];

%ignore voxels parallel to fixel of interest
idx_non_parallel=find(abs(fix_ang_sep_mat(:,4))<90-22.5);

%vox gives number of coherently oriented fixels in voxels approximately
%inline with fixel of interest

vox=sum(sum(fix_ang_sep_mat(idx_non_parallel,5:end)<seperation_between_voxels));

%if number of coherently oriented fixels less than number to compare,
%delete fixel of interest
if vox<=number_to_compare
    %delete fixel (to be done after this loop officially)
    main_fix_cart(mm+1,1:4)=[0, 0, 0, 0];
    %update number of fixels in voxel
    main_fix_cart(1,4)=main_fix_cart(1,4)-1;
end

end

%remove all rows of zeros 
main_fix_cart=main_fix_cart(any(main_fix_cart,2),:);
output_fixel=main_fix_cart;
end