function D=Dice(X,Y)

% Dice - compute the dice metric between binaries images

XnY=X .* Y;

D=2*sum(XnY(:))/(sum(X(:))+sum(Y(:)));

end

