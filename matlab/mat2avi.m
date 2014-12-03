function mat2avi(I,filename)

for i=1:size(I,4)
    imshow(I(:,:,:,i));
    f(i)=getframe;
end
%movie(f);
movie2avi(f,filename,'compression','None','quality',100);


end
