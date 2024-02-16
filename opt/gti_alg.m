function out = gti_alg(Y,prms)
   if nargin < 2
       N = size(Y,1);
       K = floor(0.1*N^2);
   else
       if isfield(prms,'N'); N = prms.N; 
       else; N = size(Y,1); end

       if isfield(prms,'K'); K = prms.K; 
       else; N = floor(0.1*N^2); end

       if isfield(prms,'gamma'); gamma = prms.gamma; 
       else; gamma = 0.1; end
   end
     
    [O,M] = size(Y);
    H = N-O;
    E = N*(N-1)/2;

    %generate L, A, B for a fully connected graph
    L = eye(N)*N-ones(N);
    A = ones(N)-eye(N);
    B = incidence(graph(A));
    
    %reconstruct Y adding zeros
    Yr = Y;%[Y;zeros(H,M)];
    %p = [ones(O,1) zeros(H,1)];
    %P = diag(p);
    In = eye(N);
        
    cvx_begin quiet
        variable w(E,1) 
        variable Z(M,M)
        expression Lw(N,N)
        expression H(N+M,N+M)
    
        
        for e = 1:E
            Lw = Lw+w(e)*B(:,e)*B(:,e)';
        end
        %Z = Yr'*inv(In+gamma*Lw)*Yr + gamma*Yr'*Lw*Y;
        H = [Z, Y'; Y,In+gamma*Lw];
    
        minimize(trace(Z))
    
        subject to
            H == semidefinite(N+M);
            ones(1,E)*w == K;
            w >= 0;
            w <= 1;
    
    cvx_end

    L_hat = zeros(N);
    [~,b] = sort(w);
    c = zeros(size(w));
    c(b(1:K)) = 1;
    for e = 1:E
            L_hat = L_hat+c(e)*B(:,e)*B(:,e)';
    end
    out.L_hat = L_hat;
end