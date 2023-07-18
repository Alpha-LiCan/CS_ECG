function [ s_OMP ] = CS_StOMP( y,Theta,epsilon )
%Orthogonal Matching Pursuit Algorithm (OMP)
%   Inputs    
%       y     : measured vector
%       Theta : Sensing matrix
%       epsilon : tolerance for approximation between successive solutions.
%   Output
%       s_OMP  : Solution found by the algorithm

%By MohammadReza Jabbari-Email: Mo.re.jabbari@gmail.com.

N = size(Theta,2);
s_OMP = zeros(N,1);
r = y;
I =[];
while norm(r,2)>=epsilon
    m = max(abs(Theta'*r));
    I = union(find(abs(Theta'*r)== m),I);
    s_OMP(I) = pinv(Theta(:,I))*y;
    r = y-Theta*s_OMP;
end 


end

