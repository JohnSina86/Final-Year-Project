clc
clear all
close all

e_height = 10.0 ; 
e_thickness = 2.0 ; 
e_width = 5.0 ; 

letter(1,1) = 0.0 ; 
letter(2,1) = e_height ; 
letter(1,2) = 0.0  ; 
letter(2,2) = 0.0 ; 
letter(1,3) = e_width ; 
letter(2,3) = 0.0 ;
letter(1,4) = e_width ; 
letter(2,4) = e_thickness ;
letter(1,5) = e_thickness ; 
letter(2,5) = e_thickness ;
letter(1,6) = e_thickness ; 
letter(2,6) = e_thickness ;
letter(1,7) = e_thickness ; 
letter(2,7) = 2*e_thickness ;
letter(1,8) = e_width ; 
letter(2,8) = 2*e_thickness ;
letter(1,9) = e_width ; 
letter(2,9) = 3*e_thickness ;
letter(1,10) = e_thickness ; 
letter(2,10) = 3*e_thickness ;
letter(1,11) = e_thickness ; 
letter(2,11) = 4*e_thickness ;
letter(1,12) = e_width ; 
letter(2,12) = 4*e_thickness ;
letter(1,13) = e_width ; 
letter(2,13) = 5*e_thickness ;
letter(1,14) = 0.0 ; 
letter(2,14) = 5*e_thickness ;
letter(3,1:14) = 1.0 ; 

theta = pi/4.0 ; 
rotation  = [ cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0 ; 0 0 1 ] ; 
translation =  [ 1 0 20 ; 0 1 -5 ; 0 0 1 ] ; 
shear  = [ 1 2 0 ; 0 1 0 ; 0 0 1 ] ; 
scale = [ 0.5 0 0 ; 0 2 0 ; 0 0 1 ] ; 
reflection = [ -1 0 0 ; 0 1 0 ; 0 0 1] ; 

transformation = reflection*scale*shear*translation*rotation ; 

new_poly = transformation*letter ; 

figure
fill(letter(1,:),letter(2,:),'k' ) ; 
hold on
fill(new_poly(1,:),new_poly(2,:),'r' ) ; 
axis([-20 20 -20 20 ]) 
