function plot_seeds(I,seeds)


[s1,s2,s3] = size(I);
colorIm = zeros(s1,s2,3,s3);

color=[1,0,0;0,1,0];

% Copy the image and boundaries into colored images space
for k=1:s3
    colorIm(:,:,:,k) = repmat( I(:,:,k), [1,1,3] );
end

Nseeds=size(seeds,2);
for label=1:2
for i=1:Nseeds
    s1=round(squeeze(seeds(label,i,:)));
    colorIm(s1(1),s1(2),:,s1(3))=color(label,:);
end
end

implay(colorIm);

%s1=squeeze(seeds(label,i,:));
%z1=round(s1(3));
%figure;
%imshow(squeeze(I(:,:,z1)));
%hold on
%plot(s1(1),s1(2),'+');
%hold off

end