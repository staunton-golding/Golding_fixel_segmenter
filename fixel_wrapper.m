function [directions_scaling, number_of_fixels, fifty_one_voxel] = fixel_wrapper(dodec_points,test_odf,model_odf,empty_mat, cap_thres,single_fiber_max)
% Usage: [directions_scaling, number_of_fixels, fifty_one_voxel] = fixel_wrapper(dodec_points,test_odf,model_odf,empty_mat, cap_thres,single_fiber_max)
% 
% fixel_wrapper segments an input odf into it's consitituent fixels. It does
% so by iteratively fitting a single-fiber FODF to the targeted 'test' FODF to a
% cylindrically symmetric set of points / cap on the test FODF, removing the
% rotated and scaled single fiber FODF each time. This function stops the
% iterative process once no cylindrically symmetric cap is able to be
% found, or once the remaining found fixels would be less than 0.1 units. This 0.1 threshold 
% standard best practice (https://onlinelibrary.wiley.com/doi/abs/10.1002/hbm.22099) 
% The default threshold of max fixel size = 0.1 units is inline with current FODF segmentation
% method standards, changing this threshold hold value is done by changing
% the last value on line 115 
% 'while abs_max>>0.1'
% 
% Input Variables:
%   
%   dodec_points: an nx3 matrix of sampled points around a unit sphere
%
%   test_odf: SH coefficient representation of FODF to be segmeneted.
%             Standard FODF SH representation (only even order FODFs.) Segments FODFs
%             of order 8.
%
%   model_odf: SH coefficient representation of single-fiber model ODF (all
%              non m=0 coefficients should = 0, to maintain cylindrical symmetry
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
%   note* model_odf, cap_thres, and single_fiber_max, are all
%   outputs of the function "single_fiber_odf" and empty_mat is output of
%   the function "sh_by_matrix"
%
% Output Variables:
%   
%   directions_scaling: nx4 matrix where columns 1:3 correspond to xyz
%                       directions of n found fixels, and column 4
%                       corresponds to scaling factor used to fit
%                       single_fiber FODF for each fixel
%
%   number_of_fixels: integer referring to number of found fixels
%
%   fifty_one_voxel: the number of times while segmenting voxel of interest
%                    a cap with 50% agreement is needed (while max fixel
%                    height is still above dual threshold)



%get ODF of interest
test_odf=squeeze(test_odf);

%get xyz coords of sampled angles
x_ang=dodec_points(:,1);
y_ang=dodec_points(:,2);
z_ang=dodec_points(:,3);

%convert cartesian coords to spherical coords (of sampled angles)
[azi_start, ela_start, ~]=cart2sph(x_ang,y_ang,z_ang);

%load model single fiber FOD 
z_sh_sym=model_odf;

%sample test-odf in equidistant directions
if size(test_odf,2)==1
    test_odf_sh=empty_mat*test_odf;
else
    test_odf_sh=empty_mat*test_odf';
end

%convert test ODF spherical values to cartesian
[x_start,y_start,z_start] = sph2cart(azi_start,ela_start,test_odf_sh);

%fixel_lib is the radius and cartesian values of the target ODF at each
%sampled direction
fixel_lib=zeros(length(azi_start),4);
fixel_lib(:,1)=test_odf_sh;
fixel_lib(:,2)=x_start;
fixel_lib(:,3)=y_start;
fixel_lib(:,4)=z_start;

%the inital azimuth and elevation angles are saved
first_sph=zeros(length(azi_start),2);
first_sph(:,1)=azi_start;
first_sph(:,2)=ela_start;

%empty arrays to save the direction of each fixel, as well as the scaling
%factor of each fixel
inds_for_fixels=zeros(1,100);
c_all=zeros(1,100);

%counter
aa=0;
fifty_one_voxel=0;
%initalize a max absolute value greater than 0.005 


%the overall maximum value of the FODF, with no fixel removed
starting_max=max(max(fixel_lib(:,1)));
abs_max=starting_max;

while abs_max>0.1
    %increase counter (for number of fixels found)
    aa=aa+1;
    
    %kernel around which to construct fixel (max radial value)
    %ind_fixel_dir is index of fixel direction used later with original ela
    %and azi 
    kern=(max(max(fixel_lib(:,1))));
    [ind_fixel_dir]=find(fixel_lib==kern);
    ind_fixel_dir=ind_fixel_dir(1);
    
    %get azimuth of max value (kernel) (elevation found after first
    %rotation)
    az=first_sph(ind_fixel_dir, 1);
 
    %vectorize x, y, z values - need to vectorize from fixel_lib bc
    %fixel_lib updated with each deleted fixel
    x_mat=fixel_lib(:,2);
    x_vec=x_mat(:);
    
    y_mat=fixel_lib(:,3);
    y_vec=y_mat(:);
    
    z_mat=fixel_lib(:,4);
    z_vec=z_mat(:);
    
    %put old xyz into matrix - for rotation
    xyzold=[x_vec, y_vec,z_vec]';
    
    %find rotation around z and y axis, using the azimuth and elevation of
    %max angle - this brings the max value to parallel with the positive z
    %axis, first rotate around z with azimuth to bring kernel into xz plane
    rz=rotz(rad2deg(-az));
     
    xyz_first_rot=rz*xyzold;
      
    %seperate rotated xyz matrix
    x_first_rot=xyz_first_rot(1,:);
    y_first_rot=xyz_first_rot(2,:);
    z_first_rot=xyz_first_rot(3,:);

    %now find the elavation needed for rotation around y axis
    [~,ela_rot1,~]=cart2sph(x_first_rot,y_first_rot,z_first_rot);
    ela_rot1_use=ela_rot1(ind_fixel_dir);

    %rotate around y axis
    ry=roty(rad2deg((pi/2+ela_rot1_use)));
    xyz_new=ry*xyz_first_rot;
    x_new=xyz_new(1,:);
    y_new=xyz_new(2,:);
    z_new=xyz_new(3,:);
    
    %establish the absolute max = kernel, and the index is the same as
    %above, abs_max also defined at end of loop after fixel removed, so
    %that a new cycle wont be started with a too small abs_max
    abs_max=kern;

    for delt_eigen=20:5:50
    %parameter for difference in delta eigenvalues (percent ratios)

    %initialize stat to be greater than delt_eigen, and start looking at cap = all
    %values > 0.5 abs max 
    stat=55;
    dr=0.5*abs_max;
    delta_shell=0.01;

    while stat>delt_eigen
        
        %step the shell up towards the max value
        dr=dr-delta_shell;

        %use z axis to evaluate shell inclusion (helps for excluding seperate radial
        %peaks, as their z axis would be low(er))
        shell_circ=find(z_new<=abs_max&z_new>abs_max-dr);
        
        %identify the new circular cap / shell
        x_shell=x_new(shell_circ);
        y_shell=y_new(shell_circ);
        z_shell=z_new(shell_circ);
        
        %construct second moment matrix of shell/cap
        M=zeros(2,2);
        M(1,1)=sum(z_shell.*(x_shell.^2));
        M(1,2)=sum(z_shell.*x_shell.*y_shell);
        M(2,1)=sum(z_shell.*x_shell.*y_shell);
        M(2,2)=sum(z_shell.*(y_shell.^2));
        
        %find eigenvalues of above matrix, evaluate to determine difference
        %in eigenvalues
        e=eig(M);
        stat=abs(e(1)-e(2))/(e(1)+e(2))/2 * 100;
    end
   
    %find spherical coordinate representation of cylindrically symmetric
    %shell/cap
    [azi_shell, ela_shell, rad_shell]=cart2sph(x_shell,y_shell,z_shell);

    %find spherical coordinate rep of rotated cartesian coordinates (whole
    %FODF)
    [azi_rot, ela_rot, radi]=cart2sph(x_new,y_new,z_new);

    %only seclect azimuth and elevation within cylindrically symmetric
    %shell/cap and within monotonic cap of single fiber ODF
    idx_shell=(abs(ela_shell)>cap_thres);
  
    ela_shell=ela_shell(idx_shell);
    azi_shell=azi_shell(idx_shell);
    rad_shell=rad_shell(idx_shell);
    
    %if shell insufficient in size, break from while loop
    if length(rad_shell)>3
       break 
    end
   
    end
    
    %if shell sufficient in size, use it
    if length(rad_shell)<3
        break
    end

    %if lowest amount of agreement (50%) used while neither max value threshold
    %reached, denote this
    if delt_eigen==50
     fifty_one_voxel=fifty_one_voxel+1;
    end

    %ind_fixel_dir is the index of peak of fixel - used for orientation evaluation
    %later
    inds_for_fixels(aa)=ind_fixel_dir;
    
    %evaluate single fiber odf both only at angles correlating to shell,
    %and at angles equal to the rotated test_odf
    
    sh_z_norm=spherical_reconstruction(z_sh_sym,8,azi_shell+pi,ela_shell+pi/2);
    sh_z_fit=spherical_reconstruction(z_sh_sym,8,azi_rot+pi,ela_rot+pi/2);

    %Initalize linear system of equations
    prob=eqnproblem;
    c=optimvar('c',1);
    opts=optimoptions(@lsqlin,'Display','off');
    for jj=1:length(ela_shell)
        %set up each specific equation (shell in direction = scaling factor times *
        %sampled single fiber ODF in that direction)
        x=string(jj);
        prob.Equations.(append('eqn',x)=rad_shell(jj)==(c(1))*sh_z_norm(jj);
    
    end
    
    %initial guess at solution
    c0.c=1;
    %solve system of equations
    [sol,~,~]=solve(prob,c0,'Options',opts);
    
    %identify scaling factor
    fac=sol.c;
    %record scaling factor for future use
    c_all(aa)=fac;

    %scale the whole single fiber odf by the correct factor
    sh_z_fit=sh_z_fit*fac;
    
    %subtract scaled single fiber FODF from target FODF, prevent negative
    %values (for any overestimation by fit)
    sh_z_new=radi-sh_z_fit;
    sh_z_new(sh_z_new<0)=0;

    %convert scaled single fiber odf to cartesian coordinates
    [x_rot_del,y_rot_del,z_rot_del] = sph2cart(azi_rot,ela_rot,sh_z_new);
  
    %find new max, and repopulate fixel_lib with updated values
    fixel_lib(:,1)=(x_rot_del.^2+y_rot_del.^2+z_rot_del.^2).^.5;
    fixel_lib(:,2)=x_rot_del;
    fixel_lib(:,3)=y_rot_del;
    fixel_lib(:,4)=z_rot_del;
    
    %reset abs_max
    abs_max=(max(max(fixel_lib(:,1))));

    %first_sph needs to be rewritten - so that new fixel removed FOD can be
    %reoriented with new max value (could use original azi and ela angle,
    %but would require rewriting new matrix of fixel-removed FOD in
    %cartesian coordinates, clunky)
    [azi_rot_del, ela_rot_del, ~]=cart2sph(x_rot_del,y_rot_del,z_rot_del);
   
    first_sph=zeros(length(azi_rot_del),2);
    first_sph(:,1)=azi_rot_del;
    first_sph(:,2)=ela_rot_del;
    
end

%remove nonzero inputs into indices for fixels
inds_for_fixels=nonzeros(inds_for_fixels);
c_all=nonzeros(c_all);

%get unit vectors for each direction (could use azi_start and ela_start,
%doing this because 3.14~=0 with cart2sph and so rounding errors can cause
%slight issues it seems)
[azi_true,ela_true,~]=cart2sph(x_start,y_start,z_start);
[x_fix, y_fix, z_fix]=sph2cart(azi_true, ela_true, ones(length(azi_true)));

%normalize cartesian coords of each fixel where each fixel has magnitutde =
%abs value of fixel
scaling=c_all.*single_fiber_max;

%put directions and scaling into output vector
directions_scaling=[x_fix(inds_for_fixels),y_fix(inds_for_fixels),z_fix(inds_for_fixels),scaling];
number_of_fixels=size(c_all,1);
end
