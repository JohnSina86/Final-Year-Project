% Code to implement some basic matrix operations 

clc
clear all
close all

 fprintf(1,'Lets look at matrix multiplication:  \n') ; 

% A matrix can be created as follows
 A = [ 2 3 1  ; 2 0 -4 ]  % This is a 2*3 matrix
 v = [1  ;4 ;5 ] % This is a 3*1 matrix (i.e. column vector)
 product = A*v   
 
  fprintf(1,'Press any key to continue \n') ; 
 pause

 % Alternatively we can use for loops to create the matrix 
 no_of_rows = 3 ; 
 no_of_cols = 5 ; 
 for(ct1 = 1:no_of_rows ) 
 for(ct2 = 1:no_of_cols ) 
 B(ct1,ct2) = rand() ;  % row ct1 col ct2 has random entry 
 end
 end
 
 fprintf(1,'Lets look at identifying sub-matrices   \n') ;

 B 
 % We can specify sub-matrices as follows
 sub_matrix_of_B = B(1:2,2:4) % Extract sub-matrix - rows 1 to =2 and columns 2 to 4  
 
 fprintf(1,'Press any key to continue \n') ; 
  pause

 fprintf(1,'Lets look at transposes \n') ; 
 
 % We can transpose matrices using the ' operator 
 B_transpose = B'
 
 fprintf(1,'Press any key to continue \n') ; 
  pause

 for(ct1 = 1:no_of_rows ) 
 for(ct2 = 1:no_of_cols ) 
 B_complex(ct1,ct2) = rand() + i*rand() ;  % row ct1 col ct2 has random entry 
 end
 end
 
 fprintf(1,'Lets look at complex matrices  \n') ; 
 
  B_complex 
  
  % Note that the ' operator produces the conjugate transpose
 B_conjugate_transpose = B_complex'

 fprintf(1,'Press any key to continue \n') ; 
  pause

  % To get the unconjugated transpose 
  B_transpose  = conj(B_complex)' 
  
   fprintf(1,'To get the inner product   \n') ; 
   a = [ 1 3  4 ] 
   c = [ 2 ; -3 ; 5] 
   inner_prod = a*c 
   new_matrix = c*a 
   
   C = rand(5,5)  
   
   the_inverse_of_C = inv(C)  
   C*the_inverse_of_C 
   the_inverse_of_C*C
   
   simple_identity = eye(5) 
