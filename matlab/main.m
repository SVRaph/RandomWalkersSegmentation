%% Segmentation by retrieval with guided random walks: 
%  Application to left ventricle segmentation in MRI
%  d'après A. Eslami, A. Karamalis, A. Katouzian, N. Navab

%% Rozand, Sarrazin, Sivera
%  2014.11.18

%% Hypothèses
%   - 3D images
%   - drivers are in a nx(sz) matrix
%   - patient image is of size sz, 
%   - myocardia have the same center and their size are normalized
%   - same mean in each image

%% Paramètres
Dice_threshold=0.5;
seg_threshold=0.5; % Threshold to turn the result into binaries
Nseeds = 6;

alpha = 90 /100;
beta  = 130/100;
gamma = 0.4;


%% Lecture des données
sz = [100 100 100]/2;
c = [50 50 30]/2;
ax1 = [1 1 0];
ax2 = [-1 1 0];
ax3 = [0 0 1];
axes = [ax1/sqrt(sum(ax1.*ax1)) ; ax2/sqrt(sum(ax2.*ax2)) ; ax3/sqrt(sum(ax3.*ax3))];
rays = [40 20 5]/2;

[R,B,~] = generate_heart(sz, c, axes, rays,0.3,Nseeds,3);

R=reshape(R,[1,sz]);
B=reshape(B,[1,sz]);

[I0,B0,seeds] = generate_heart(sz, c, axes, rays,0.3,Nseeds,3);


%% Initialisation
D_max=-1;
indx=-1;
X_opt=zeros(sz);


%% ROI -- NOT USED
for i=[]
if 0
    g_fg=squeeze(mean(seeds(1,:,:)));
    g_bg=squeeze(mean(seeds(2,:,:)));

    dx=128;%max(x)-min(x);
    dy=128;%max(y)-min(y);

    %roi=min(512,max(1,[xg-dx,xg+dx,yg-dy,yg+dy]));
    roi=[xg-dx,xg+dx,yg-dy,yg+dy];
    if (any(roi<1) || any(roi>512))
        fprintf('Your argument is invalid\n');
        return;
    end

    I=I0(roi(1):roi(2),roi(3):roi(4));
    I=I+moyR-mean(I);
end
end

I=I0;

%% Main loop on drivers
for k=1:size(R,1)
    R_k=squeeze(R(k,:,:,:));
    B_k=squeeze(B(k,:,:,:));
    X_k=Guided_Random_Walks(I,R_k,B_k,seeds,alpha,beta,gamma);
    D_k=Dice(X_k,B_k);
    
    if D_k>D_max
        D_max=D_k;
        indx=k;
        Xopt=X_k;
    end
    
end
 
if D_max>Dice_threshold
    fprintf(['Best segmentation with driver ',int2str(indx),'\nDice metrix=',num2str(D_max),'\n'])
    %implay(Xopt);
else
    fprintf('No matching subject found\nPerforming conventional Random Walks\n');
    %Xopt=Random_Walks(I);  
end

% Display
implay(B_k,5);
implay((Xopt>seg_threshold),5);

color_result = show_boundaries(I0,B0,Xopt,seg_threshold);
