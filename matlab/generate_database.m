%% Generate data base

sz       = [100 100 100]; % taille des images
n_drivers= 10;            % nombre de référence
sizeblur = 3;

% coeur moyen
center = [50 50 50];
ax1 = [1 1 0];
ax2 = [-1 1 0];
ax3 = [0 0 1];
axis = [ax1/sqrt(sum(ax1.*ax1)) ; ax2/sqrt(sum(ax2.*ax2)) ; ax3/sqrt(sum(ax3.*ax3))];
rays = [15 25 10];
width= 0.3;

%s = RandStream('mcg16807','Seed',0);

% génération des coeurs
%[R,B,~]=generate_heart( sz, center, axis, rays, width, 0,sizeblur);
%implay(R);


R=zeros([n_drivers,sz]);
B=zeros([n_drivers,sz]);
for k=1:n_drivers
    centerk= center+0.05*rand(1,3).*sz;
    axisk  = axis + 0.1*rand(3,3);
    raysk  = rays + 0.05*rand(1,3).*sz;
    
    [Rk,Bk,~]=generate_heart( sz, centerk, axisk, raysk, width, 0,sizeblur);
    
    R(k,:,:,:)=Rk;
    B(k,:,:,:)=Bk;
    disp(k);
end


% génération du patient
[I,B0,seeds]=generate_heart( sz, center, axis, rays, width, 12,sizeblur);
%implay(R);

save('data.mat','R','B','I','B0','seeds');
