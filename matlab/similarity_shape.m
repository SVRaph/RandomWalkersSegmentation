function S = similarity_shape(I,R,B,indx,X_opt,seg_threshold)


%% Création de la distance à la segmentation de référence 
SpeedImage = ones([100 100 100]);
SourcePoint=[];

R_k=squeeze(R(indx,:,:,:));
B_k=squeeze(B(indx,:,:,:));
X_k = (X_opt > seg_threshold)*1;

for i=1:size(B_k,1);
    for j=1:size(B_k,2);
        for k=1:size(B_k,3);
            if B_k(i,j,k)==1
                SourcePoint(:,end+1)=[i;j;k];
            end
        end
    end
end
LX = msfm3d(SpeedImage,SourcePoint,true,true);

%% Calcul des gradients des images Raw et test
sigma=1;
[Rx,Ry,Rz]=gradient(R_k);
Rnorm=sqrt(Rx.^2+Ry.^2+Rz.^2);

[Ix,Iy,Iz]=gradient(I);
Inorm=sqrt(Ix.^2+Iy.^2+Iz.^2);

%% Calcul du coefficient de similarité de forme
dX = 0;
S = 0;
for i=1:size(R_k,1);
    for j=1:size(R_k,2);
        for k=1:size(R_k,3);
            if X_k(i,j,k)==1
                if X_k(i+1,j,k)*X_k(i-1,j,k)*X_k(i,j+1,k)*X_k(i,j-1,k)*X_k(i,j,k+1)*X_k(i,j,k-1)==0
                    dX = dX+1;
                    scal=Rx(i,j,k)*Ix(i,j,k)+Ry(i,j,k)*Iy(i,j,k)+Rz(i,j,k)*Iz(i,j,k);
                    S_upd = abs(scal)/(Rnorm(i,j,k)+eps(1))/(Inorm(i,j,k)+eps(1))*exp(-LX(i,j,k)^2/sigma);
                    S = S + S_upd;
                end
            end
        end
    end
end

fprintf('Shape similarity :\n');
S=S/dX

end

