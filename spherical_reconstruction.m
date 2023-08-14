function sh_r=spherical_reconstruction(slm, order, azi, ela)
% Usage: sh_r=spherical_reconstruction(slm, order, azi, ela)
% 
% spherical_reconstruction_normal2 evaluates the spherial harmonic
% expansion given a set of coefficients at the inputted azimuth/elevation
% pairs. Although it is a relatively fast function, it is far quicker to
% multiply the coefficients times the previously written 'empty_mat' (look
% at sh_by_matrix).
%
% Input Variables:
%
%   slm: even ordered spherical harmonic coefficients 
%
%   degree: max degree to which the spherical harmonic will be evaluated
%
%   azi: the azimuth angles
%
%   ela: the inclination angles (harmonicY uses inclination --> 0 to pi
%   measured from z axis, not elevation --> -pi/2 to pi/2 from +x axis)
%
% Output Variables:
%
%   sh_r: The sum of all evaluated spherical harmonics times their
%         respective coefficients at the sampled angle pairs
%
% note* calls function 'harmonicY'

        pp=0;
    for kk=0:order/2
        l=kk*2;
        
        %loop through all possible orders
        for m=-l:l
            %increase coefficient counter
            pp=pp+1;
       
            
            % compute REAL spherical harmonic of degree l and order m
            if m<0
                
                Y = (2^.5)*slm(pp)*imag(harmonicY(l,(-1*m),ela,azi));
                
            elseif m==0
              
                Y= slm(pp)*harmonicY(l,m,ela,azi);
                
            else
                
                Y= (2^.5)*slm(pp)*real(harmonicY(l,(m),ela,azi));
               
            end
            
            %add SH cumulatively through 45 coefficients (addedfun
            %is the FODF)
            if pp==1
                addedfun=Y;
            else
                addedfun=addedfun+Y;
            end 
        end
 
    end

   sh_r=addedfun; 
end