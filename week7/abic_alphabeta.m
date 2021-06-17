function abic = abic_alphabeta(d,W,G,alpha,L,beta)
% function to calculate A-BIC value given 
% d - data 
% G - design matrix for a linear problem
% W - weighting function (generally inverse of data covariance matrix)
% alpha - smoothing hyper-parameter
% L - smoothing matrix (generally finite-difference Laplacian operator) 
% beta - length hyper-parameter
% OUTPUT
% abic - information criteria
% Rishav Mallick, EOS, 2019
% reference Fukuda and Johnson (2008; 2010), Funning et al., (2014)


npatch = length(L);
m = (G'*W*G + (alpha^2)*(L'*L) + (beta^2)*eye(npatch))\(G'*W*d);

abic1 = length(d)*log((d-G*m)'*W*(d-G*m) + (alpha^2)*(L*m)'*(L*m) + (beta^2)*(m'*m));

eigprior = abs(eig((alpha^2)*(L'*L) + (beta^2)*eye(npatch) ));
eigprior(eigprior<=0) = [];
abic2 = sum(log(eigprior));

eigdes = abs(eig(G'*W*G + (alpha^2)*(L'*L) + (beta^2)*eye(npatch) ));
eigdes(eigdes<=0) = [];
abic3 = sum(log(eigdes));

% return ABIC value
abic = abic1 - abic2 + abic3;

end