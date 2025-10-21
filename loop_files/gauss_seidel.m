% Code to apply Gauss Seidel to a matrix equation obtained by applying the Method
% of Moments to a wave scattering problem

clc
clear all
close all

% Get the matrix from another code 
TX = -5 ;
TY = 10 ;  
length_of_strip  = 10.0 ; 
the_data = make_the_mom_matrix(length_of_strip,TX,TY) ;

temp = size(the_data)  ;   
N = temp(1) ;  
A = the_data(1:N,1:N) ; 
b = the_data(1:N,N+1)  ;
the_locations = the_data(1:N,N+2) ; 

exact_x = inv(A)*b ; 

no_of_iterations = 2 ; 
x_GS = zeros(N,1) ;
  
for(iteration_counter = 1:no_of_iterations) 
    for(ct1 = 1:N)  
        sum = b(ct1)  ;
        for(ct2 = 1:ct1-1)
          sum = sum - A(ct1,ct2)*x_GS(ct2); 
        end
        for(ct2 = ct1+1:N)
          sum = sum  - A(ct1,ct2)*x_GS(ct2) ; 
        end
        x_GS(ct1) = sum / A(ct1,ct1)  ; 
    end
end

figure 
plot(the_locations,abs(x_GS)) 
hold on
plot(the_locations,abs(exact_x),'r-.') 
grid on
legend('GS','Exact') ; 
title(['Current on plate after ' num2str(no_of_iterations) ' iterations of GS'])
xlabel('distance on plate')
ylabel('abs( current)') 
