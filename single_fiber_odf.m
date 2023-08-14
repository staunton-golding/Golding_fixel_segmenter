function [model_odf,cap_thres,single_fiber_max]=single_fiber_odf(wmfod_matlab,voxels_matlab,points,empty_mat)
% Usage: [model_odf,cap_thres,single_fiber_max]=single_fiber_odf(wmfod_matlab,voxels_matlab,points,empty_mat)
%
% single_fiber_odf uses the FODFs from the voxels used in WM response
% function estimation to approximate an ODF of a uniformly oriented fiber
% population (a single fixel). FODFs from the aforementioned voxels are
% oriented along the z axis and their rotated coefficients are averaged.
% All non m=0 coefficients are made = 0 to preserve cylindrical symmetry
%
% note* this function calls the function 'getSHrotMtx' 
%
% Input variables:
% 
%    wmfod_matlab: 4D matrix where the first three dimensions correspond to xyz
%                 indices and the 4th dimension corresponds to SH coeffecients
%                 of FODFs (45 coefficients per FODF)
%   voxels_matlab: 4D matrix where the first three dimensions correspond to xyz
%                  indices and the 4th dimension is essentially a 3 layered
%                  mask, where the third input corresponds to whether that
%                  voxel was used in response function estimation of white
%                  matter (first two inputs correspond to rf of gm and csf)
%
%   points: an nx3 matrix of sampled points around a unit sphere, where
%           the columns are xyz
% 
%   empty_mat: nx45 matrix where index (n,m) correlates to the nth spherical
%              harmonic expansion at point m, where each column is an
%              evaluated spherical harmonic, stored in 'standard' order:
%              (l,m) = (0,0), (1,-1), (1,0), (1,1), (2,-1), ...
%
% Output Variables:
%
%   model_odf: SH coefficient representation of single-fiber model ODF (all
%              non m=0 coefficients should = 0, to maintain cylindrical symmetry
%
%   cap_thres: The elevation wherein, if starting at max radial value, the
%              single-fiber model FODF loses its monotonicity
%
%   single_fiber_max: the max radial value of the single-fiber model ODF
%
% note* calls the function 'getSHrotMtx' 

%get xyz values of equidsitant sampled points
x_ang=points(:,1);
y_ang=points(:,2);
z_ang=points(:,3);

%convert points to spherical representation
[azi_start, ela_start, ~]=cart2sph(x_ang,y_ang,z_ang);

%voxels used in wm_rf calculation
wm_data=voxels_matlab(:,:,:);

%indices of voxels used in wm_rf calculation
[indx, indy, indz]= ind2sub(size(wm_data),find(wm_data ~= 0));

%preallocate array for SH coefficients of single fiber model ODF
z_sh=zeros(45,1);

for jj=1:length(indx)
%evaluate spherical harmonic representation of FODF of voxel used in wm_rf
%calculation
single_fiber_odf=wmfod_matlab(indx(jj),indy(jj),indz(jj),:);
single_fiber_odf=squeeze(single_fiber_odf);
if size(single_fiber_odf,2)==1
    sh_single_fiber=empty_mat*single_fiber_odf;
else
    sh_single_fiber=empty_mat*single_fiber_odf';
end

%find max value of FODF
kern=(max(max(sh_single_fiber)));

%find index of max value of FODF
[ind1]=find(sh_single_fiber==kern);
ind1=ind1(1);


%find elevation and azimuth of max radius value 
az=azi_start(ind1);
el=ela_start(ind1);

%rotation around z and y axis
rz=rotz(rad2deg(-az));
ry=roty(rad2deg(pi/2-el));

%rotation matrix
rotmat=ry*rz;

%spherical harmonic rotation matrix
rotmat_sh=getSHrotMtx(rotmat,8,'real');

%fill in odd valued SH coefficients (all = 0) and create logical vector of
%even/odd coefficients
sh_zero=single_fiber_odf(1)';
sh_zero_logic=ones(1);

sh_one=zeros(1,3);

sh_two=single_fiber_odf(2:6)';
sh_two_logic=ones(1,length(sh_two));

sh_three=zeros(1,7);

sh_four=single_fiber_odf(7:15)';
sh_four_logic=ones(1,length(sh_four));

sh_five=zeros(1,11);

sh_six=single_fiber_odf(16:28)';
sh_six_logic=ones(1,length(sh_six));

sh_seven=zeros(1,15);

sh_eight=single_fiber_odf(29:45)';
sh_eight_logic=ones(1,length(sh_eight));

sh_logical_vector=[sh_zero_logic,sh_one,sh_two_logic,sh_three,sh_four_logic,sh_five,sh_six_logic,sh_seven,sh_eight_logic]';
sh_original_vec=[sh_zero,sh_one,sh_two,sh_three,sh_four,sh_five,sh_six,sh_seven,sh_eight]';

%get spherical harmonic representation of FODF used in wm_rf calculation
%oriented parallel to z axis
z_norm_sh=rotmat_sh*sh_original_vec;

%keep only even degree SH coefficients
z_norm_sh=z_norm_sh(sh_logical_vector==1);

%normalize coefficients to have max coefficient = 1
z_norm_sh=z_norm_sh./(max(z_norm_sh));

%add coefficients (later to be divided as to get average odf)
z_sh=z_sh+z_norm_sh;

end

%average SH coefficients
z_sh=z_sh./jj;

%only use all non m=0 entries (so that model ODF is cylindrically symmetric
%about z axis)
model_odf=zeros(1,45);
model_odf(1)=z_sh(1);
model_odf(4)=z_sh(4);
model_odf(11)=z_sh(11);
model_odf(22)=z_sh(22);
model_odf(37)=z_sh(37);

%identify elevation where single fiber ODF loses monotonicity, and max of
%single fiber ODF radially
azi_cap=zeros(1,500);
ela_cap=linspace(-pi/2,pi/2,500);

capper=spherical_reconstruction_normal2(model_odf,8,azi_cap,ela_cap+pi/2);

[x,~,z] = sph2cart(azi_cap,ela_cap,capper);
single_fiber_max=max(max(z));
cap_thres=abs(ela_cap(x==max(x)));

end