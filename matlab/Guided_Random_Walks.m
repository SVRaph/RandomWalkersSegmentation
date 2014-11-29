function X_k = Guided_Random_Walks(I,R,B,seeds,alpha,beta,gamma)


% Guided_Random_Walk - segmente l'image I en utilisant le driver (R,Bcl)
sz=size(I);
v_sub2ind=[1,sz(1),sz(1)*sz(2)];
N=sz(1)*sz(2)*sz(3);

disp(['Starting guided random walks N=',int2str(N)]);

% 6-voisinage
sub_c_II= [1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1]              + 2;
c_II =    sub2ind(sz,sub_c_II(:,1),sub_c_II(:,2),sub_c_II(:,3)) - sub2ind(sz,2,2,2) ; % pdt scalR with v ?
c_IR = [c_II;0]; %ind_c_IR=[c_II;0,0,0];

    function [W,O,L,A]=weights_matrices(I,R,N,c_II,c_IR,alpha,beta)
        % compute the W and the Omega matrices
        % W=spalloc(N,N,size(c_II,1)*N);
        % O=spalloc(N,N,size(c_IR,1)*N);
        
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
        
        % Omega
        nzerosmax=size(c_IR,1)*N;
        v=zeros(nzerosmax,3);
        vsum=zeros(N,1);
        k=1;
        for i=1:N
            for l=1:size(c_IR,1)
                j=i+c_IR(l);
                if (j >= 1 && j <= N)
                    v(k,:)=[i,j,exp(-beta*(I(i)-R(j))^2)];
                    vsum(i)=vsum(i)+v(k,3);
                    k=k+1;
                end
            end
        end
        O=sparse(v(1:k-1,1),v(1:k-1,2),v(1:k-1,3),N,N,nzerosmax);
        A=spdiags(vsum,0,N,N); % Ai=sum(O(i,:)) 
        
    end



% Matrices W, Omega, L, A
disp('Compute W, Omega, L and A');
[W,O,L,A]=weights_matrices(I,R,N,c_II,c_IR,alpha,beta);
disp('L,A');
disp('done');

% Indices Marked et Unmarked
indM1=reshape(seeds(1,:,:)-1,size(seeds,2),3)*v_sub2ind'+1;
indM2=reshape(seeds(2,:,:)-1,size(seeds,2),3)*v_sub2ind'+1;
indM=[indM1;indM2];

logicalM=false(N,1);
for i=1:size(indM,1)
    logicalM(round(indM(i)))=true;
end
logicalU=not(logicalM);

% Système sparse
Lu=L(logicalU,logicalU);
Lb=L(logicalM,logicalU);
Au=A(logicalU,logicalU);
Ab=A(logicalM,logicalU);
Omegab=O(logicalM,logicalU);
Omegau=O(logicalU,logicalU);

bm=B(logicalM);
bu=B(logicalU);
xm=[zeros(size(indM1,1),1);ones(size(indM1,1),1)];


MA = Lu+gamma*Au;
Mb = (gamma/2)*Omegab'*bm + (gamma/2)*Omegau*bu - Lb'*xm - gamma*Ab'*xm;

disp('Sparse linear system created');

% solve MA*x=Mb
xu=pcg(MA,Mb);

X_k=zeros(size(B));
X_k(logicalU)=xu;

disp('Sparse linear system solved');

end