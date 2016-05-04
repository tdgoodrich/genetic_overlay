function [L,D] = build_Laplacian(W)

D = diag(sum(W));

L = D - W;

end









