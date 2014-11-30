function X_k = Random_Walks(I,seeds,alpha)


% Guided_Random_Walk - segmente l'image I en utilisant le driver (R,Bcl)
sz=size(I);
v_sub2ind=[1,sz(1),sz(1)*sz(2)];
N=sz(1)*sz(2)*sz(3);

disp(['Starting guided random walks N=',int2str(N)]);

% 6-voisinage
sub_c_II= [1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1]              + 2;
c_II =    sub2ind(sz,sub_c_II(:,1),sub_c_II(:,2),sub_c_II(:,3)) - sub2ind(sz,2,2,2) ; % pdt scalR with v ?


% 26-voisinage
%sub_c_II= [1,0,0 ; -1,0,0 ; 0,1,0 ; 1,1,0 ; -1,1,0 ; 0,-1,0 ; 1,-1,0 ; -1,-1,0;...
%           1,0,1 ; -1,0,1 ; 0,1,1 ; 1,1,1 ; -1,1,1 ; 0,-1,1 ; 1,-1,1 ; -1,-1,1; 0,0,1 ...
%           1,0,-1 ; -1,0,-1 ; 0,1,-1 ; 1,1,-1 ; -1,1,-1 ; 0,-1,-1 ; 1,-1,-1 ; -1,-1,-1; 0,0,-1] +2 ;
%c_II =    sub2ind(sz,sub_c_II(:,1),sub_c_II(:,2),sub_c_II(:,3)) - sub2ind(sz,2,2,2) ; % pdt scalR with v ?

    function [W,L]=weights_matrices(I,N,c_II,alpha)
        % compute the W and the Omega matrices
        
        % W
        nzerosmax=size(c_II,1)*N;
        v=zeros(nzerosmax,3);
        vsum=zeros(N,1);
        k=1;
        for i=1:N
            for l=1:size(c_II,1)
                j=i+c_II(l);
                if (j >= 1 && j <= N)
                    v(k,:)=[i,j,exp(-alpha*(I(i)-I(j))^2)];
                    vsum(i)=vsum(i)+v(k,3);
                    k=k+1;
                end
            end
        end
        W=sparse(v(1:k-1,1),v(1:k-1,2),v(1:k-1,3),N,N,nzerosmax);
        L=-2*W+spdiags(2*vsum,0,N,N);      
    end



% Matrices W, Omega, L, A
fprintf('Compute W, L...');
[W,L]=weights_matrices(I,N,c_II,alpha);
fprintf(' done \n');

% Indices Marked et Unmarked
indM1=reshape(seeds(1,:,:)-1,size(seeds,2),3)*v_sub2ind'+1;
indM2=reshape(seeds(2,:,:)-1,size(seeds,2),3)*v_sub2ind'+1;
indM=round([indM1;indM2]);

logicalM=false(N,1);
for ii=1:size(indM,1)
    logicalM(indM(ii))=true;
end
logicalU=not(logicalM);

% Système sparse
Lu=L(logicalU,logicalU);
Lb=L(logicalU,logicalM);

xm=[zeros(size(indM1,1),1);ones(size(indM2,1),1)];


MA = Lu;
Mb = - Lb*xm;

disp('Sparse linear system created');

% solve MA*x=Mb
%xu=pcg(MA,Mb,1e-6,100);
xu=pcg(MA,Mb,1e-5,500);

X_k=zeros(size(I));
X_k(logicalU)=xu;
X_k(logicalM)=xm;
disp('Sparse linear system solved');

end