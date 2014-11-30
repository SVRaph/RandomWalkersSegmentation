function [ colorIm ] = show_boundaries (I0,B0,Bk)

% show the image with the theoretical and computed segmentation boundaries
[s1,s2,s3] = size(I0);
colorIm = zeros(s1,s2,3,s3);

% Compute the theoretical boundaries
dif1 = B0-circshift(B0,[1,1]);
dif2 = B0-circshift(B0,[1,2]);
notTheoBound = (abs(dif1)+abs(dif2)==0);

% Compute the computed boundaries
dif1 = Bk-circshift(Bk,[1,1]);
dif2 = Bk-circshift(Bk,[1,2]);
notCompBound = (abs(dif1)+abs(dif2)==0);

% Copy the image and boundaries into colored images space
for k=1:s3
    theoBound = ones(s1,s2)-notTheoBound(:,:,k);
    compBound = ones(s1,s2)-notCompBound(:,:,k);
    colorIm(:,:,:,k) = repmat( I0(:,:,k) .* notTheoBound(:,:,k) .* notCompBound(:,:,k) , [1,1,3])...
        + cat(3 , theoBound.*notCompBound(:,:,k) , zeros(s1,s2) , zeros(s1,s2))...
        + cat(3 , zeros(s1,s2) , compBound , zeros(s1,s2));
end

implay(colorIm);

end
