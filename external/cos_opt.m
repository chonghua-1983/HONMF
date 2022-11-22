function W = cos_opt(X)
% input
% X:term-document matrix, row represents document

W = X*X';
DX = sqrt(diag(W))*sqrt(diag(W))';
W = W./DX;

end
