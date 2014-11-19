function [ I ] = gaussian_blur_3d(I,n)

% flou gaussien isotrope sur une image 3d
% sigma=n/7
if n==0
    return;
end

sigma=n/7;
x = ( (0:n-1)-(n-1)/2 );
f = exp( -x.^2/(2*sigma^2) );
f = f / sum(f(:));

Gx=zeros(n,1,1);
Gy=zeros(1,n,1);
Gz=zeros(1,1,n);

Gx(:,1,1)=f;
Gy(1,:,1)=f;
Gz(1,1,:)=f;

I = convn(I,Gx,'same'); %convolution in x
I = convn(I,Gy,'same'); %convolution in y
I = convn(I,Gz,'same'); %convolution in z
end