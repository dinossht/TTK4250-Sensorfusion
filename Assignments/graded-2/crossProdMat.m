function S = crossProdMat(n)
% It takes a vector n ∈ R3 and creates the skew symetric matrix that 
% corresponds to the matrix multiplication implementation of the cross 
% product. 
% Book: equation (10.5)

    S = [0      -n(3)   n(2)    ;
         n(3)   0       -n(1)   ;
         -n(2)  n(1)    0       ;];
end