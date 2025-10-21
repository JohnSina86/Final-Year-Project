% Code to illustrate the application of affine transformations to 2D polygons
clc
clear all
close all

%   Define a 3*5 matrix made up of five 3*1 column vectors 
%   These vectors are the vertices of the polygon. 

polygon_point(1,1) = 0.0 ;  
polygon_point(2,1) = 0.0 ; 
polygon_point(3,1) = 1.0 ; 
polygon_point(1,2) = 10.0  ; 
polygon_point(2,2) = 0.0 ;  
polygon_point(3,2) = 1.0 ;
polygon_point(1,3) = 10.0 ; 
polygon_point(2,3) = 10.0 ; 
polygon_point(3,3) = 1.0 ;
polygon_point(1,4) = 0.0;  
polygon_point(2,4) = 10.0 ; 
polygon_point(3,4) = 1.0 ;
polygon_point(1,5) = polygon_point(1,1) ;  
polygon_point(2,5) = polygon_point(2,1) ;  
polygon_point(3,5) = 1.0 ;

centre_x = 0.25*(polygon_point(1,1) + polygon_point(1,2) + polygon_point(1,3) + polygon_point(1,4)) ;
centre_y = 0.25*(polygon_point(2,1) + polygon_point(2,2) + polygon_point(2,3) + polygon_point(2,4)) ;

theta = 3*pi/4.0  ; 

rot_matrix  = [ cos(theta)  -sin(theta) 0 ; sin(theta) cos(theta)  0 ; 0 0 1 ]  ;
rot_polygon  = rot_matrix*polygon_point ; 

tx = 2.0 ; 
ty = 4.0 ; 

trans_matrix  = [ 1  0 tx  ;  0 1 ty ;  0 0 1  ]  ;
trans_polygon  = trans_matrix*polygon_point ; 

ref_matrix  = [ -1  0 0  ;  0 1 0 ;  0 0 1  ]  ;
ref_polygon  = ref_matrix*polygon_point ; 

scaling_matrix  = [ 0.5  0 0  ;  0 0.5 0 ;  0 0 1  ]  ;
scaled_polygon  = scaling_matrix*polygon_point ; 

shear_matrix  = [ 1  0.5 0  ;  0 1 0 ;  0 0 1  ]  ;
shear_polygon  = shear_matrix*polygon_point ; 

scaled_then_translated_polygon = trans_matrix*scaling_matrix*polygon_point ; 

figure
plot(polygon_point(1,1:5),polygon_point(2,1:5),'-*') 
grid on
hold on
fill(rot_polygon(1,1:5),rot_polygon(2,1:5),'g-*')  
plot(trans_polygon(1,1:5),trans_polygon(2,1:5),'y-*')  
plot(ref_polygon(1,1:5),ref_polygon(2,1:5),'k-*')  
plot(scaled_polygon(1,1:5),scaled_polygon(2,1:5),'m-*')  
plot(shear_polygon(1,1:5),shear_polygon(2,1:5),'r-*')  
plot(scaled_then_translated_polygon(1,1:5),scaled_then_translated_polygon(2,1:5),'c-*')  
axis([-30 30 -30 30]) 
legend('Original','Rotated','Translated','Reflected','scaled','shear','scaled then translated') ; 
