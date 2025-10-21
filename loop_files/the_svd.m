% A code to implement the SVD 
clc 
clear all
close all

%   Specify the number of rows and columns 
N = 6 ; 
M = 8 ; 

% Make a random N*M matrix
% Could use A = rand(N,M) 
% Prefer to explicitly make the matrix 

for(ct1 = 1:N) 
   for(ct2  = 1:M)
       A(ct1,ct2) = rand() ; 
   end
end

A 
[U,S,V] = svd(A) 
recreated_A = U*S*V'
