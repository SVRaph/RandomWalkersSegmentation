function [maxDice, maxIndex,abscisse] = test_parameters(I,B0,R,B,seeds,alpha, beta, gamma, threshold,drivers)

sz=size(R);
sz=sz(2:end);

sa = size(alpha,2);
sb = size(beta,2);
sg = size(gamma,2);
sv = max( max(sa,sb),sg) ;

maxDice = zeros( sv,1 );
maxIndex= zeros( sv,1 );
n=1;

abscisse = ((sa>1)*alpha + (sb>1)*beta + (sg>1)*gamma)';

for a=1:sa
    for b=1:sb
        for g=1:sg

%% Initialisation
D_max=-1;
indx=-1;
X_opt=zeros(sz);

%% Main loop on drivers
for k=drivers
    R_k=squeeze(R(k,:,:,:));
    B_k=squeeze(B(k,:,:,:));
    X_k=Guided_Random_Walks(I,R_k,B_k,seeds,alpha(a),beta(b),gamma(g));
    D_k=Dice(X_k>0.36,B_k);
    
    if D_k>D_max
        D_max=D_k;
        indx=k;
        X_opt=X_k;
    end
    fprintf(['Dice index ',num2str(D_k),'\n']);
end

maxDice(n) = Dice(X_opt>threshold, B0)
maxIndex(n)= indx;
n=n+1;
        end
    end
end


end