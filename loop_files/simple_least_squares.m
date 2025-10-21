%   Code to implement simple least squares 
clc 
clear all 
close all

% Read in the data values 
data_points = load('data.res')
x = data_points(:,1) ; 
y =  data_points(:,2) ;
the_N = size(x) 
N = the_N(1)  

sum_xy = 0.0  ; 
sum_xx = 0.0 ; 

for(ct = 1:N ) 
sum_xy = sum_xy + x(ct)*y(ct) ; 
sum_xx = sum_xx + x(ct)*x(ct) ;
end

b_tilde = sum_xy / sum_xx ; 
fprintf(1,'b_tilde is %1.10f \n',b_tilde) ; 

for(ct = 1:N) 
A(ct,1) = x(ct) ; 
y_vec(ct,1) = y(ct) ; 
end

answer =  inv(A'*A)*A'*y_vec ; 
alt_b  = answer ; 
fprintf(1,'Using Linear Algebra we get %1.10f \n',alt_b) ; 

not_as_good_b = b_tilde*0.99  ; 

error = 0.0 ;
not_as_good_error = 0.0 ; 
for(ct = 1:N) 
alt_predicted_y(ct) =  alt_b*x(ct) ; 
predicted_y(ct) = b_tilde*x(ct) ; 
not_as_good_y(ct) = not_as_good_b*x(ct) ; 
error = error + power(y(ct) - predicted_y(ct),2.0) ;
not_as_good_error = not_as_good_error + power(y(ct) - not_as_good_y(ct),2.0) ;
end

fprintf(1,'Min error equals %f \n',error) ; 
fprintf(1,'Not as good error equals %f \n',not_as_good_error) ; 

figure
plot(x,y,'*')
hold on
plot(x,predicted_y,'r') ;
plot(x,alt_predicted_y,'g-.') ;
plot(x,not_as_good_y,'k') ;
title('Least Squares Example') 
xlabel('x')
ylabel('y') 
grid on 
legend('data','LS1','LS2','not as good') ;
