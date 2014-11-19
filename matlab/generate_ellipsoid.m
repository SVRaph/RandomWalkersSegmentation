function [ e ] = generate_ellipsoid( s, c, X, R)
% generate a matrix in M_{s[1] x s[2] x s[3]} containing a quadratic function
% for which each iso-surface is an ellipsoid with parameters:
%   - c vector in R^3     : center of the ellipsoid
%   - X matrix in R^(3x3) : directions of the ellipsoid
%   - R vector in R^3     : rays of the ellipsoid

% How to generate parameters: (CTRL+T the next line and change the values
% s = [s1 s2 s3];
% c = [c1 c2 c3];
% X1 = [x11 x12 x13];
% X2 = [x21 x22 x23];
% X3 = [x31 x32 x33];
% X = [X1/sqrt(sum(X1.*X1)) ; X2/sqrt(sum(X2.*X2)) ; X3/sqrt(sum(X3.*X3))];
% R = [R1 R2 R3];

% generate whole picture with zeros
e = zeros(s(1), s(2), s(3));

% generate probabilities whether the voxel is inside or outside the
% ellipsoid

Sigma= X'*diag(R.^(-2))*X;

for i=1:s(1)
    for j=1:s(2)
        for k=1:s(3)
            voxel = ([i j k]-c)';
            e(i,j,k) = voxel' * Sigma * voxel;
        end
    end
end

end

