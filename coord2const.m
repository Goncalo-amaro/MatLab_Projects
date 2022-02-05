function [C] = coord2const( X, w )
%COORD2CONST - Calculation of the HCW eq. constants 
%   

C(1) = 2*X(3)+X(4)/w;
C(2) = X(6)/w;
C(3) = -3*X(3)-2*X(4)/w;
C(4) = X(1)-2*X(6)/w;
C(5) = X(5)/w;
C(6) = X(2);

end

