function [ R,B,seeds ] = generate_heart( sz, center, axis, rays, width, Nseeds,sizeblur)

% generate a synthetic heart, its segmentation and some seeds

% USAGE: [R_k,B_k,seeds_k] = generate_heart(s, c, X, R,0.3,3); implay(R_k);

% intput 
%   - sz     in N^3     : size of the image
%   - center in R^3     : center of the ellipsoid
%   - axis   in R^(3x3) : directions of the ellipsoid
%   - rays   in R^3     : rays of the ellipsoid
%   - width  in R       : parametrise the width of the heart muscle
%   - Nseeds in N       : number of seeds to generate

% output
%   - R 3D image of size sz
%   - B 3D binary image of size sz (segmentation)
%   - seeds in N^(2xNseedsx3) : seeds(1,k,:) voxel coordinates of the k_th background seed  


% Images
e = generate_ellipsoid(sz,center,axis,rays);

B = (abs(e-1)<width);

noise=0.1*randn(size(B));
bias =zeros(size(B));

R = gaussian_blur_3d(B,sizeblur)+ noise + bias; % blur(B)


% Seeds
seeds=zeros(2,Nseeds,3);

M =  axis'*diag(rays);
for k=1:2:Nseeds
    r=randn(3,1);
    r=r/norm(r);
    seeds(1,k,:)=center'+sqrt(1-2*width)*M*r;
end
for k=2:2:Nseeds
    r=randn(3,1);
    r=r/norm(r);
    seeds(1,k,:)=center'+sqrt(1+2*width)*M*r;
end
for k=1:2:Nseeds
    r=randn(3,1);
    r=r/norm(r);
    seeds(2,k,:)=center'+sqrt(1-0.4*width)*M*r;
end
for k=2:2:Nseeds
    r=randn(3,1);
    r=r/norm(r);
    seeds(2,k,:)=center'+sqrt(1+0.4*width)*M*r;
end



end


