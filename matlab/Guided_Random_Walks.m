function X_k = Guided_Random_Walks(I,R,B,seeds,alpha,beta,gamma)

%% TODO 
% - commentaires
% - correction erreurs présentes dans l'article
% - efficaité

% Guided_Random_Walk - segmente l'image I en utilisant le driver (R_k,B_k)

sz=size(I);
v_sub2ind=[1,sz(1),sz(1)*sz(2)];
N=sz(1)*sz(2)*sz(3);

disp(['Starting guided random walks ',int2str(N)]);

% 6-voisinage
sub_c_II= [1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1]              + 2;
c_II =    sub2ind(sz,sub_c_II(:,1),sub_c_II(:,2),sub_c_II(:,3)) - sub2ind(sz,2,2,2) ; % pdt scalR with v ?
c_IR = [c_II;0]; %ind_c_IR=[c_II;0,0,0];

    function [W,O]=weights_matrices(I,R,N,c_II,c_IR,alpha,beta)
        % compute the W and the Omega matrices
        W=spalloc(N,N,size(c_II,1)*N);
        O=spalloc(N,N,size(c_IR,1)*N);
        for i=1:N
            if not(mod(i,100))
                disp(i);
            end
            for l=1:size(c_II,1)
                j=i+c_II(l);
                if (j >= 1 && j <= N)
                    W(i,j) = exp(-alpha*(I(i)-I(j)));
                end
            end
            
            for l=1:size(c_IR,1)
                j=i+c_IR(l);
                if (j >= 1 && j <= N)
                    O(i,j) = exp(-beta*(I(i)-R(j)));
                end
            end
        end
    end

    function [L,A]=energy_matrices(W,O)
        L=-2*W;
        A=spalloc(N,N,N);
        
        for i=1:N
            if not(mod(i,100))
                disp(i);
            end
            L(i,i)=2*sum(W(i,:));
            A(i,i)=sum(O(i,:));
        end
    end


% Matrices W, Omega, L, A
disp('W');
[W,O]=weights_matrices(I,R,N,c_II,c_IR,alpha,beta);
disp('L,A');
[L,A]=energy_matrices(W,O);
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
Mb = gamma*Omegab'*bm + gamma*Omegau*bu - Lb'*xm + gamma*Ab'*xm;

disp('Sparse linear system created\n');

% solve MA*x=Mb
xu=pcg(MA,Mb);

X_k=zeros(size(B));
X_k(logicalU)=xu;

disp('Sparse linear system solved\n');

end