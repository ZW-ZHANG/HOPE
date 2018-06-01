function [UL, thetas, VL, WL] = jdgsvds( A, B, L, tol, K)
[N N1] = size(A);
[N2 N3] = size(B);
if N1 ~= N || N2 ~= N || N3 ~= N
    fprintf('Error: the matrix must be square\n');
end
if nargin < 5
    K = L;
end
U = zeros(N, K);
W = zeros(N, K);
c = zeros(K,1);
d = zeros(K,1);
r = zeros(2*N,1);

% GMRES variables
inK = 30;
Q = zeros(2*N,inK+1); 
tmpQ = zeros(2*N,1);
Hes = zeros(inK+1, inK);

% sigular vectors and values
% L = 1;
thetas = zeros(L);
UL = zeros(N, L);
WL = zeros(N, L);

% buffers
AW = zeros(N, K);
ATU = zeros(N, K);
BBW = zeros(N, K);
H = zeros(K, K);
ATUL = zeros(N, L);
BBWL = zeros(N, L);
s = rand(N,1);
t = rand(N,1);
k = 0;
theta = 0.0;
for l = 1:L
    if k < 1
        s = s./sqrt(s'*s);
        bt = B*laptrans(WL(:,1:l), BBWL(:,1:l)', t);
        t = t./sqrt(bt'*bt);
        U(:,1) = s;
        W(:,1) = t;
        k = 1;
        theta = 0.0;
    end
%     k = 1;
%     for k = 1:2
%         AW(:, k) = laptrans(UL(:,1:l), UL(:,1:l)', A*laptrans(WL(:,1:l), BBWL(:,1:l)', W(:,k)));
%         ATU(:, k) = laptrans(BBWL(:,1:l), WL(:, 1:l)', A'*laptrans(UL(:,1:l), UL(:,1:l)', U(:,k)));
%         BBW(:, k) = laptrans(BBWL(:,1:l), WL(:,1:l)', B'*(B*laptrans(WL(:,1:l), BBWL(:,1:l)', W(:,k))));
%         H(k,1:k) = U(:,k)'*AW(:,1:k);
%         H(1:k,k) = ATU(:,1:k)'*W(:,k);
%         
%         t = laptrans(BBWL(:,1:l), WL(:,1:l)', A'*laptrans(UL(:,1:l), UL(:,1:l)', laptrans(UL(:,1:l), UL(:,1:l)', A*laptrans(WL(:,1:l), BBWL(:,1:l)', t))));
%         W(:,k+1) = t - W(:,1:k)*(t'*BBW(:,1:k))';
%         bw = B*(W(:,k+1)-WL(:,1:l)*(BBWL(:,1:l)'*W(:,k+1)));
%         W(:,k+1) = W(:,k+1)./sqrt(bw'*bw);
%         U(:,k+1) = laptrans(UL(:,1:l), UL(:,1:l)', laptrans(UL(:,1:l), UL(:,1:l)', A*laptrans(WL(:,1:l), BBWL(:,1:l)', W(:,k+1))));
%         U(:,k+1) = laptrans(U(:,1:k), U(:,1:k)', U(:,k+1));
%         U(:,k+1) = U(:,k+1)./sqrt(U(:,k+1)'*U(:,k+1));
%     end
    bg = clock();
    while true
        AW(:, k) = laptrans(UL(:,1:l), UL(:,1:l)', A*laptrans(WL(:,1:l), BBWL(:,1:l)', W(:,k)));
        ATU(:, k) = laptrans(BBWL(:,1:l), WL(:, 1:l)', A'*laptrans(UL(:,1:l), UL(:,1:l)', U(:,k)));
        BBW(:, k) = laptrans(BBWL(:,1:l), WL(:,1:l)', B'*(B*laptrans(WL(:,1:l), BBWL(:,1:l)', W(:,k))));
        H(k,1:k) = U(:,k)'*AW(:,1:k);
        H(1:k,k) = ATU(:,1:k)'*W(:,k);
    %     H = U(:,1:k)'*A*W(:,1:k);
        pre_theta = theta;
        [c(1:k,1), theta, d(1:k,1), tmpflag] = svds(H(1:k,1:k),1);
        u = U(:,1:k)*c(1:k,1);
        w = W(:,1:k)*d(1:k,1);
        bbw = BBW(:,1:k)*d(1:k,1);
        r(1:N) = AW(:,1:k)*d(1:k,1)-theta*u; 
        r((N+1):(2*N)) = ATU(:,1:k)*c(1:k,1)-theta*bbw;
        ed = clock();
        %fprintf('%d-th singular value, %d-th iteration, theta: %f, error: %f %f %d, time: %fs\n', l, k, theta, sqrt(r(1:N)'*r(1:N)), sqrt(r'*r), tmpflag, etime(ed,bg));
        bg = clock();
        if sqrt(r'*r) <= tol || K <= 1
            UL(:,l) = u;
            WL(:,l) = w;
            BBWL(:,l) = B'*(B*WL(:,l));
            ATUL(:,l) = A'*UL(:,l);
            thetas(l,l) = theta;
            
            resn = min(max([10,floor(2*K/3)]), k);
            if resn > 1
                [resc, resthe, resd] = svds(H(1:k,1:k), resn);
                U(:,1:(resn-1)) = U(:,1:k)*resc(:,2:resn);
                W(:,1:(resn-1)) = W(:,1:k)*resd(:,2:resn);
                AW(:,1:(resn-1)) = AW(:,1:k)*resd(:,2:resn);
                ATU(:,1:(resn-1)) = ATU(:,1:k)*resc(:,2:resn);
                BBW(:,1:(resn-1)) = BBW(:,1:k)*resd(:,2:resn);
                H(1:(resn-1),1:(resn-1)) = U(:,1:(resn-1))'*AW(:,1:(resn-1));
                theta = H(1,1);
            end
            k = resn-1;
            break;
        end

        if k == K
            resn = max([min([K-5,10]), floor(2*K/3)]); 
            [resc, resthe, resd] = svds(H(1:K,1:K), resn);
            U(:,1:resn) = U(:,1:K)*resc;
            W(:,1:resn) = W(:,1:K)*resd;
            AW(:,1:resn) = AW(:,1:K)*resd;
            ATU(:,1:resn) = ATU(:,1:K)*resc;
            BBW(:,1:resn) = BBW(:,1:K)*resd;
            H(1:resn,1:resn) = U(:,1:resn)'*AW(:,1:resn);
            k = resn;
        end

        if (abs(pre_theta-theta) < 0.001*abs(pre_theta+theta))
            %GMRES
            beta = sqrt(r'*r);
            Q(:,1) = -r./beta;
            rgmres = 1.0;
            for ink = 1:inK
%                 tmpQ(1:N) = A*(Q((N+1):(2*N),ink)-WL(:,1:l)*(BBWL(:,1:l)'*Q((N+1):(2*N),ink)));
%                 tmpQ(1:N) = tmpQ(1:N)-UL(:,1:l)*(UL(:,1:l)'*tmpQ(1:N));
%                 tmpQ(1:N) = tmpQ(1:N)-theta*Q(1:N,ink);
% %                 tmpQ((N+1):(2*N)) = A'*(Q(1:N,ink)-UL(:,1:l)*(UL(:,1:l)'*Q(1:N,ink)))-theta*(B'*(B*(Q((N+1):(2*N),ink)-WL(:,1:l)*(BBWL(:,1:l)'*Q((N+1):(2*N),ink)))));
%                 tmpQ((N+1):(2*N)) = A'*Q(1:N,ink)-ATUL(:,1:l)*(UL(:,1:l)'*Q(1:N,ink))-theta*(B'*(B*Q((N+1):(2*N),ink))-BBWL(:,1:l)*(BBWL(:,1:l)'*Q((N+1):(2*N),ink)));
%                 tmpQ((N+1):(2*N)) = tmpQ((N+1):(2*N))-BBWL(:,1:l)*(WL(:,1:l)'*tmpQ((N+1):(2*N)));
                
                tmpQ(1:N) = A*Q((N+1):(2*N),ink)-theta*Q(1:N,ink);
%                 tmpQ((N+1):(2*N)) = A'*Q(1:N,ink)-theta*B'*(B*Q((N+1):(2*N),ink));
                tmpQ((N+1):(2*N)) = (Q(1:N,ink)'*A)'-theta*((B*Q((N+1):(2*N),ink))'*B)';

%                 tmpQ(1:N) = tmpQ(1:N,1)-u*(u'*tmpQ(1:N,1))-UL*(UL'*tmpQ(1:N,1));
                tmpQ(1:N) = tmpQ(1:N,1)-u*(u'*tmpQ(1:N,1))-UL*(tmpQ(1:N,1)'*UL)';
%                 tmpQ((N+1):(2*N)) = tmpQ((N+1):(2*N),1)-bbw*(w'*tmpQ((N+1):(2*N),1))-BBWL*(WL'*tmpQ((N+1):(2*N),1));
                tmpQ((N+1):(2*N)) = tmpQ((N+1):(2*N),1)-bbw*(w'*tmpQ((N+1):(2*N),1))-BBWL*(tmpQ((N+1):(2*N),1)'*WL)';
%                 Hes(1:ink,ink) = Q(:, 1:ink)'*tmpQ;
                Hes(1:ink,ink) = (tmpQ'*Q(:, 1:ink))';
                Q(:, ink+1) = tmpQ-Q(:,1:ink)*Hes(1:ink,ink);
                Hes(ink+1,ink) = sqrt(Q(:, ink+1)'*Q(:, ink+1));
                Q(:, ink+1) = Q(:,ink+1)./Hes(ink+1,ink);
                pre_rgmres = rgmres;
                rgmresv = Hes(1:(ink+1),1:ink)*((Hes(1:(ink+1),1:ink)'*Hes(1:(ink+1),1:ink))\(beta*Hes(1,1:ink)'))-beta*eye(ink+1,1);
                rgmres = sqrt(rgmresv'*rgmresv);
                if  abs(pre_rgmres-rgmres) < 0.0001*(pre_rgmres+rgmres)
                    break;
                end
            end
            tmpQ = Q(:,1:ink)*((Hes(1:(ink+1),1:ink)'*Hes(1:(ink+1),1:ink))\(beta*Hes(1,1:ink)'));
            s = tmpQ(1:N,1);
            t = tmpQ((N+1):(2*N),1);
        elseif abs(pre_theta-theta) < 0.1*abs(pre_theta+theta) && l==1 && false
            Minv = [ones(N,1)./(theta+0.1); sum(B.^2,1)'./(theta+0.1)];
            tmpm = [(u.*Minv(1:N))'*u 0; 0 (w.*Minv((N+1):(2*N)))'*bbw];
            Minvr = Minv.*r;
            tmpst = [Minv(1:N).*u zeros(N,1); zeros(N,1) Minv((N+1):(2*N)).*bbw]*(tmpm\[u'*Minvr(1:N);w'*Minvr((N+1):(2*N))]);
            tmpst = Minvr-tmpst;
            s = tmpst(1:N);
            t = tmpst((N+1):(2*N));
        else
            s = r(1:N);
            t = r((N+1):(2*N));
        end
        s = laptrans(UL(:,1:l), UL(:,1:l)', s);
        t = laptrans(WL(:,1:l), BBWL(:,1:l)', t);

        U(:,k+1) = s - U(:,1:k)*(s'*U(:,1:k))';
        U(:,k+1) = U(:,k+1)./sqrt(U(:,k+1)'*U(:,k+1));
        W(:,k+1) = t - W(:,1:k)*(t'*BBW(:,1:k))';
        bw = B*(W(:,k+1)-WL(:,1:l)*(BBWL(:,1:l)'*W(:,k+1)));
        W(:,k+1) = W(:,k+1)./sqrt(bw'*bw);
        k = k+1;
    end
    
    if k < 1
        s = rand(N,1);
        t = rand(N,1);
        s = laptrans(UL(:,1:l), UL(:,1:l)', s);
        t = laptrans(WL(:,1:l), BBWL(:,1:l)', t);
    end
end
VL = B*WL(:,1:L);
    function rVec = laptrans(M1, M2, Vec)
        rVec = Vec - M1*(M2*Vec);
    end
    function y = correction_multiply(x)
        y = zeros(2*N,1);
        y(1:N) = A*(x((N+1):(2*N))-WL(:,1:l)*(BBWL(:,1:l)'*x((N+1):(2*N))));
        y(1:N) = y(1:N)-UL(:,1:l)*(UL(:,1:l)'*y(1:N));
        y(1:N) = y(1:N)-theta*x(1:N);
%                 y((N+1):(2*N)) = A'*(Q(1:N,ink)-UL(:,1:l)*(UL(:,1:l)'*Q(1:N,ink)))-theta*(B'*(B*(Q((N+1):(2*N),ink)-WL(:,1:l)*(BBWL(:,1:l)'*Q((N+1):(2*N),ink)))));
        y((N+1):(2*N)) = A'*x(1:N)-ATUL(:,1:l)*(UL(:,1:l)'*x(1:N))-theta*(B'*(B*x((N+1):(2*N)))-BBWL(:,1:l)*(BBWL(:,1:l)'*x((N+1):(2*N))));
        y((N+1):(2*N)) = y((N+1):(2*N))-BBWL(:,1:l)*(WL(:,1:l)'*y((N+1):(2*N)));

        y(1:N) = y(1:N,1)-u*(u'*y(1:N,1));
        y((N+1):(2*N)) = y((N+1):(2*N),1)-bbw*(w'*y((N+1):(2*N),1));
    end
end

