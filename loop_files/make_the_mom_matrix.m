function result = make_the_mom_matrix(length_of_strip,TX,TY) 

% Don't change anything in this code! 

% Define the frequency, wavenemuber and free space parameters
MHz = 1000000.0 ; 
f = 2000.0*MHz ; 
% Assume free space (vaccuum )
epsilon0 = 8.854e-12 ; 
mu0 = 4.0*pi*1.0e-7 ; 
% Speed of light in vaccum
c0 = 1.0/sqrt(epsilon0*mu0) ; 
omega = 2.0*pi*f ; 

%Wavenumber 
k0 = omega/c0  ; 
% Wavelength
lambda = c0/f ;
%Impedance 
eta = sqrt(mu0/epsilon0) ; 

no_of_strips = 1 ; 
strip_length(1)  = length_of_strip ; 
angle(1) = 0.0 ;

strip_start_x(1) = 0.0*lambda ; 
strip_start_y(1) = 0.0 ; 

contour_plot_x_start  = -2.0 ;
contour_plot_x_finish  = 2.0 ;
contour_plot_y_start  = -2.0 ; 
contour_plot_y_finish  = 2.0 ; 
do_contour_plots  = 1 ; 

disc_per_lambda = 10 ; 
delta_s = lambda/disc_per_lambda ; 
delta_s_contour = delta_s ; 
contour_N = ceil((contour_plot_x_finish -contour_plot_x_start)/delta_s_contour ) ; 

pos_counter = 0  ; 
for(ct_strip = 1:no_of_strips) 
    N_per_strip = ceil(strip_length(ct_strip)/delta_s) ;
    if(ct_strip >1 ) 
    strip_start_x(ct_strip) = strip_start_x(ct_strip-1) + strip_length(ct_strip-1)*cos(angle(ct_strip-1)) ;   
    strip_start_y(ct_strip) = strip_start_y(ct_strip-1) + strip_length(ct_strip-1)*sin(angle(ct_strip-1)) ;   
    end
    for(ct_within_strip = 1:N_per_strip) 
        pos_counter = pos_counter + 1 ; 
        x(pos_counter) =  strip_start_x(ct_strip) + (ct_within_strip-0.5)*cos(angle(ct_strip))*delta_s ;
        y(pos_counter)  = strip_start_y(ct_strip) + (ct_within_strip-0.5)*sin(angle(ct_strip))*delta_s;
        arc_length(pos_counter) = (pos_counter-0.5)*delta_s ;
    end
end

for(ct = 1:pos_counter) 
    Rx = x(ct) - TX ; 
    Ry = y(ct) - TY ; 
    R = sqrt(Rx*Rx + Ry*Ry) ; 
    V(ct,1) = besselj(0,k0*R) - j*bessely(0,k0*R) ; 
end

self = -1.0*delta_s*(1.0 - j*(2.0/pi)*log(1.781*k0*delta_s/(4.0*exp(1.0)))) ; 
for(ct1 = 1:pos_counter) 
    for(ct2 = 1:pos_counter) 
        Rx = x(ct1)  - x(ct2) ; 
        Ry = y(ct1)  - y(ct2) ; 
        R = sqrt(Rx*Rx + Ry*Ry) ;
        the_hank = -(besselj(0,k0*R) - j*bessely(0,k0*R)) ; 
        if( ct1 ~= ct2) 
            Z(ct1,ct2) = (k0*eta/4.0)*the_hank*delta_s ;
        else
            Z(ct1,ct2) = (k0*eta/4.0)*self ; 
        end
    end
end
    
result = [Z  V  x'] ;
