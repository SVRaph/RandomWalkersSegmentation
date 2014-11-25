function X_k = Guided_Random_Walks(I,R_k,B_k)

% Guided_Random_Walk - segmente l'image I en utilisant le driver (R_k,B_k)
% X_k vector in sz^3



X_k=B_k;

%parameters
beta = 1;
alpha = 1;
gamma = 1;

% TODO
im_size = size(I,1)*size(I,2)*Size(I,3);
vecI = reshape(I,im_size,1,1);
vecR = reshape(R_k,im_size,1,1);
vecB = reshape(B_k,im_size,1,1);

% Type of connection intra-image
connecI6 = [[0;1;0],[1;0;0],[0;0;1]];
connecI26 = [connec6,[1;1;0],[1;-1;0],...
    [1;1;1],[1;-1;1],[-1;-1;1],[-1;1;1],[1;0;1],[0;1;1],[-1;0;1],[0;-1;1];
% 26-voisinage
%sub_c_II= [1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1;...
%           1,1,0;1,-1,0;1,1,1;1,-1,1;-1,-1,1;-1,1,1;1,0,1;0,1,1;-1,0,1;0,-1,1];
%c_II =    sub2ind(sz,sub_c_II(:,1),sub_c_II(:,2),sub_c_II(:,3)) - sub2ind(sz,2,2,2) ; % pdt scalR with v ?
%c_IR = [c_II;0];

function coordN = Mat2Vec(coordV)
    coordN = coordV(1) * [1,size(I,1),size(I,2)];
end

connecR7 = [connecI6,[0;0;0]];
connecR27 = [connecI26,[0;0;0];

connecI = connecI6;
connecR = connecR7;

% compute omega
omega = zeros(im_size,im_size);
A = zeros(im_size,im_size);

for i=1:im_size
    for j=1:size(connecI,2)
        omega(i,i+Mat2Vec(connecR)) = exp(-alpha*(vecI(i)-vecR(i+Mat2Vec(connecR)))^2);
    end
end

omega = omega + omega';



end