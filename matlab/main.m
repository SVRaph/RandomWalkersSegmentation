%% Segmentation by retrieval with guided random walks: 
%  Application to left ventricle segmentation in MRI
%  d'après A. Eslami, A. Karamalis, A. Katouzian, N. Navab

%% Rozand, Sarrazin, Sivera
%  2014.11.18

%% Hypothèses
%   - 3D images
%   - drivers are in a nx(sz) matrix
%   - patient image is of size sz, 
%   - myocardia have the same center and their sizes are normalized
%   - same mean in each image

%% Paramètres
Dice_threshold=0.5; 
seg_threshold=0.36; % Threshold to turn the result into binaries

alpha = 15;
beta  = 25;
gamma = 0.4;


%% Lecture des données :
% drivers : R, B,
% patient à segmenter : I, B0 et seeds
load('data.mat'); 
sz=size(R);
sz=sz(2:end);

[~,~,seeds]=generate_heart( sz, center, axis, rays, width, 200,sizeblur);

%% Test parametres
[maxDice, maxIndex, abscisse] = test_parameters(I,B0,R,B,seeds,10, 30, 0.4, 0.36,1:size(R,1))


[maxDice, maxIndex, abscisse] = test_parameters(I,B0,R,B,seeds,10,[0.1 0.3 1 3 10 20 25 30 40 50 100 300 1000], 0.4, 0.4,7)


[maxDice, maxIndex, abscisse] = test_parameters(I,B0,R,B,seeds,15, 20, 10, 0.48,7)

semilogx(abscisse,maxDice,'-o');
xlabel('beta');
ylabel('Dice');
set(gca,'XTick',abscisse);


vseeds=[2 4 8 16 60 212];
maxDice=zeros(size(vseeds));
n=1;
for ns=vseeds
    [d, ~, ~] = test_parameters(I,B0,R,B,seeds(:,1:ns,:),15, 20, 0.4, 0.4,7);
    maxDice(n)=d;
    n=n+1;
end

semilogx(vseeds,maxDice,'-o');
xlabel('Nombre de seeds');
ylabel('Dice');
set(gca,'XTick',2*vseeds);


%% Initialisation
D_max=-1;
indx=-1;
X_opt=zeros(sz);

%% ROI -- NOT DONE

%% Main loop on drivers
for k=7%1:size(R,1)
    R_k=squeeze(R(k,:,:,:));
    B_k=squeeze(B(k,:,:,:));
    X_k=Guided_Random_Walks(I,R_k,B_k,seeds,alpha,beta,gamma);
    D_k=Dice(X_k,B_k);
    
    if D_k>D_max
        D_max=D_k;
        indx=k;
        X_opt=X_k;
    end
    fprintf(['Dice index ',num2str(D_k),'\n']);
    
end
 
if D_max>Dice_threshold
    fprintf(['Best segmentation with driver ',int2str(indx),'\nDice metrix=',num2str(D_max),'\n'])
else
    fprintf('No matching subject found\nPerforming conventional Random Walks\n');
    X_rw=Random_Walks(I,seeds,alpha);  
end

% Test against the true segmentation
% X_star=Guided_Random_Walks(I,I,B0,seeds,alpha,beta,gamma);
