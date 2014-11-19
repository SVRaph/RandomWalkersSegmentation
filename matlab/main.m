%% Segmentation by retrieval with guided random walks: 
%  Application to left ventricle segmentation in MRI
%  d'après A. Eslami, A. Karamalis, A. Katouzian, N. Navab

%% Rozand, Sarrazin, Sivera
%  2014.11.18

%% Hypothèses
%   - 2D images
%   - every driver has size 256x256 and the same mean
%   - patient image is 512x512

%% Paramètres
Dice_threshold=0.5;


%% Lecture des données -- TODO
R=zeros(256,256,5);
B=zeros(256,256,5);

moyR=mean(R(:));

I0=zeros(512,512);


%% Initialisation
D_max=-1;
indx=-1;
X_opt=zeros(size(B(:,:,1)));


%% Seeds - TODO foreground/background
imshow(I0);
[xs,ys] = ginput;

%% ROI
xg=mean(xs);
yg=mean(ys);

dx=128;%max(x)-min(x);
dy=128;%max(y)-min(y);

%roi=min(512,max(1,[xg-dx,xg+dx,yg-dy,yg+dy]));
roi=[xg-dx,xg+dx,yg-dy,yg+dy];
if (any(roi<1) || any(roi>512))
    mprintf('Your argument is invalid\n');
    return;
end

I=I0(roi(1):roi(2),roi(3):roi(4));
I=I+moyR-mean(I);


%% Main loop on drivers
for k=1:size(R,3)
    R_k=R(:,:,k);
    B_k=B(:,:,k);
    X_k=Guided_Random_Walks(I,R_k,B_k);
    D_k=Dice(X_k,B_k);
    
    if D_k>D_max
        D_max=D_k;
        indx=k;
        Xopt=X_k;
    end
    
end
    
if D_max>Dice_threshold
    mprintf('Best segmentation with driver '+int2str(indx)+'\nDice metrix='+num2str(D_max)+'\n')
    imshow(Xopt);
else
    mprintf('No matching subject found\nPerforming conventional Random Walks\n');
    Xopt=Random_Walks(I);  
end
imshow(Xopt);
