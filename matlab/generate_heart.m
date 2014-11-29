function [ R,Bg,seeds ] = generate_heart( sz, center, axis, rays, width, Nseeds,sizeblur)

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
%   - Bg 3D binary image of size sz (segmentation du ventricule gauche)
%   - seeds in N^(2xNseedsx3) : seeds(1,k,:) voxel coordinates of the k_th background seed  


% ventricule gauche
eg = generate_ellipsoid(sz,center,axis,rays);
Bg = (abs(eg-1)<width);

% ventricule droit
centerd=center+1.5*rays(1)*axis(1,:);
raysd  =1.4*rays;
ed = generate_ellipsoid(sz,centerd,axis,raysd);
Bd = (abs(eg)>(1+width)) .* (abs(ed-1)<0.5*width);

% Image bruitée
B = Bg+Bd;
noise=0.1*randn(sz);

bias =zeros(sz);
[x,y,z]=meshgrid(1:sz(1),1:sz(2),1:sz(3));
bias=bias+0.5*cos(2*x/sz(1)+3*y/sz(2)+z/sz(3));

for i=1:3
    e=generate_ellipsoid(sz,rand(1,3).*sz,rand(3,3),0.1*rand(1,3).*sz);
    bias=bias+0.4*(1+rand())*(e<1);
end

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


