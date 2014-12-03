function S = similarity_shape(B0,B_opt)


%% Création de la distance à la segmentation de référence 
SpeedImage = ones([100 100 100]);
SourcePoint=[];
B=B0;
for i=1:size(B,1);
    for j=1:size(B,2);
        for k=1:size(B,3);
            if B(i,j,k)==1
                SourcePoint(:,end+1)=[i;j;k];
            end
        end
    end
end
LB = msfm3d(SpeedImage,SourcePoint,true,true);

%% Création de la distance à la segmentation optimale
SourcePoint=[];
B=B_opt;
for i=1:size(B,1);
    for j=1:size(B,2);
        for k=1:size(B,3);
            if B(i,j,k)==1
                SourcePoint(:,end+1)=[i;j;k];
            end
        end
    end
end
LX = msfm3d(SpeedImage,SourcePoint,true,true);



%% Calcul des gradients des images LB et LX
[Rx,Ry,Rz]=gradient(LB);
Rnorm=sqrt(Rx.^2+Ry.^2+Rz.^2);

[Ix,Iy,Iz]=gradient(LX);
Inorm=sqrt(Ix.^2+Iy.^2+Iz.^2);

%% Calcul du coefficient de similarité de forme
sigma=1;
dX = 0;
S = 0;
for i=1:size(LB,1);
    for j=1:size(LB,2);
        for k=1:size(LB,3);
            if B0(i,j,k)==1
                if B0(i+1,j,k)*B0(i-1,j,k)*B0(i,j+1,k)*B0(i,j-1,k)*B0(i,j,k+1)*B0(i,j,k-1)==0
                    dX = dX+1;
                    scal=Rx(i,j,k)*Ix(i,j,k)+Ry(i,j,k)*Iy(i,j,k)+Rz(i,j,k)*Iz(i,j,k);
                    S_upd = abs(scal)/(Rnorm(i,j,k)+eps(1))/(Inorm(i,j,k)+eps(1))*exp(-LX(i,j,k)^2/sigma);
                    S = S + S_upd;
                end
            end
        end
    end
end
S=S/dX;
fprintf(['Shape similarity: ',num2str(S),'\n']);

end

