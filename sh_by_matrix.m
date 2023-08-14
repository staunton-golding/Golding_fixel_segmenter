function empty_mat=sh_by_matrix(points)
% Usage: empty_mat=sh_by_matrix(points)
% 
% sh_by_matrix evaluates all even ordered spherical harmonics up to order l=8 and their
% full ranges of degrees. The output of this function is useful in
% efficient evaluation of functions represented as spherical harmonics. 
%
% 
% Input variables:
%
%    points: an nx3 matrix of sampled points around a unit sphere, where
%           the columns are xyz. All future spherical harmonic evaluations 
%           will be limited to these points
%
% Output variables:
%
%   empty_mat: nx45 matrix where index (n,m) correlates to the nth spherical
%              harmonic expansion at point m, where each column is an
%              evaluated spherical harmonic, stored in 'standard' order:
%              (l,m) = (0,0), (1,-1), (1,0), (1,1), (2,-1), ...
%
% note* calls the function 'harmonicY'

% initialize 
empty_mat=zeros(size(points,1),45);

x_ang=points(:,1);
y_ang=points(:,2);
z_ang=points(:,3);

%convert points to spherical representationm
[azi_start, ela_start, ~]=cart2sph(x_ang,y_ang,z_ang);

%ela start is to be used with harmonicY (hence 90 deg shift up)
ela_start=ela_start+pi/2;
%azi_start is to be used with harmonicY (range needs 0-->2pi not -pi to pi)
azi_start=azi_start+pi;

%for harmonicY
degree=8;
slm=ones(1,45);

%counter through each degree/order combination (for degree 8 all 45)
pp=0;
    
    %loop through degrees
    for kk=0:degree/2
        l=kk*2;
        
        %loop through all possible orders
        for m=-l:l

            %increase coefficient counter
            pp=pp+1;
       
            % compute spherical harmonic of degree l and order m
            % (definitions same as one used in MRTrix3)

            if m<0
                
                Y = (2^.5)*slm(pp)*imag(harmonicY(l,(-1*m),ela_start,azi_start));
                
            elseif m==0
                
                Y = slm(pp)*harmonicY(l,m,ela_start,azi_start);
                
            else
                
                Y = (2^.5)*slm(pp)*real(harmonicY(l,(m),ela_start,azi_start));
               
            end
            
           %fill in each column with SH evaluation of degree kk order m
           empty_mat(:,pp)=Y;
        
    
        end
    end
    
end