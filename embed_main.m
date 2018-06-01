function [U,V] = embed_main(A, K, beta)
% Input: 
% A: N*N adjacency matrix (sparse)
% K: dimensionality of embedding space
% beta: decaying constant, default is 0.5 / spectral radius
% Output:
% U: N*K left embedding matrix
% V: N*K right embedding matrix
% The high-order proximity (katz) matrix is approximated by U * V'

[N, ~] = size(A);
% Katz: S = sum_{l=1}^{+inf}{beta*A}^l
if nargin < 3
    beta = 0.5 / getRadius(A);
end
A = beta .* A;
M = speye(N)-A;
[V, S, U] = jdgsvds(A', M', K, 0.0001, 100); % the 0.0001 error tolerance can be modified to speed up while reducing accuracy
U = U(:,1:K) * sqrt(S(1:K,1:K));
V = V(:,1:K) * sqrt(S(1:K,1:K));
end

