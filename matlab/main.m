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
Nseeds = 12;

alpha = 90 /100;
beta  = 130/100;
gamma = 0.4;


%% Lecture des données :
% drivers : R, B,
% patient à segmenter : I, B0 et seeds
load('data.mat'); 
sz=size(R);
sz=sz(2:end);


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


%% Main loop on drivers
for k=1:size(R,1)
    R_k=squeeze(R(k,:,:,:));
    B_k=squeeze(B(k,:,:,:));
    X_k=Guided_Random_Walks(I,R_k,B_k,seeds,alpha,beta,gamma);
    D_k=Dice((X_k>0.5),B_k);
    
    if D_k>D_max
        D_max=D_k;
        indx=k;
        Xopt=X_k;
    end
    fprintf(['Dice index ',num2str(D_k),'\n']);
    
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
