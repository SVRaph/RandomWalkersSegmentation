function thresh=find_seg_threshold(X_opt)

% calcul un seuil à partir d'une heuristique sur le volume

x=0.3:0.01:0.6;
vol=zeros(1,size(x,2));

for i=1:size(x,2)
    vol(i)=sum(X_opt(:)>x(i));
end
dvol=abs(vol(2:end)-vol(1:end-1));

m=mean(dvol(1:10));
i=find(dvol>3*m,1,'first');

thresh=x(i+1);

end