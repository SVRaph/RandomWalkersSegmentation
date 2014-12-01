%% Validation d'une segmentation
% Xopt against B0

seg_threshold=0.43;

find_seg_threshold(X_opt)

B1=(X_opt>seg_threshold);

% Display
%implay(abs(B1-B0));

show_boundaries(I,B0,B1);
show_boundaries(I,B0,squeeze(B(indx,:,:,:)));

% volume 
v0=sum(B0(:));
v1=sum(B1(:));
fprintf(['Relative volume error: ',num2str((v1-v0)/v0),'\n']);

% Dice
d01=Dice(B0,B1);
fprintf(['Dice index ',num2str(d01),'\n']);


%% Dice vs seg_threshold
x=0.3:0.01:0.6;
y=zeros(1,size(x,2));
vol=zeros(1,size(x,2));

for i=1:size(x,2)
    vol(i)=sum(X_opt(:)>x(i));
    y(i)=Dice(B0,X_opt>x(i));
end
%figure;plot(x,y);
%figure;plot(x,vol);

dvol=abs(vol(2:end)-vol(1:end-1));
dvol=dvol/max(dvol(:));

figure;
hold on
h=plot(x(2:end),dvol,'r'); % schéma à gauche
set(h, 'LineWidth', 2);
h=plot(x,y,'b');
set(h, 'LineWidth', 2);
hold off
xlabel('seuil de la segmentation');
legend('Dérivée normalisée du volume','Dice');


plot(x,y,'b');
xlabel('seuil de la segmentation');
ylabel('Dice');
axis('tight');

[~,i]=max(dvol)
[~,i]=max(y)




% Test against the true segmentation
X_k=Guided_Random_Walks(I,I,B0,seeds,alpha,beta,gamma);
B1=(X_k>seg_threshold);
%color_result = show_boundaries(I,B0,B1);
fprintf(['Dice index ',num2str(Dice(B0,B1)),'\n']);

% Test against the true segmentation
X_k=Guided_Random_Walks(I,I,B0,seeds,alpha,beta,0);
B1=(X_k>seg_threshold);
fprintf(['Dice index ',num2str(Dice(B0,B1)),'\n']);
