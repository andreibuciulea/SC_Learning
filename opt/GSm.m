function [L,L0] = GSm(C,reg)
    N = size(C,1);
    alpha = reg.alpha;
    beta = reg.beta;
    lambda = reg.lambda;
    Bout = gen_B12(N);
    B01 = Bout.B1;
    E = N*(N-1)/2;
    
    %disp(['alpha:' num2str(alpha) '  beta:' num2str(alpha) '  lambda:' num2str(alpha)])
    cvx_begin quiet
        variables w0(E,1)
        expressions L(N,N) symmetric semidefinite
        L = B01*diag(w0)*B01';
        minimize(alpha*vec(C)'*vec(L) + beta*vec(L(~eye(N)))'*vec(L(~eye(N)))...
            - lambda*ones(1,N)*log(diag(L)))
        subject to
            sum(abs(L*ones(N,1))) <= 1e-8;
            L(~eye(N)) <= 0;
            w0 >= 0;
    cvx_end
    
     w0 = w0>1e-5;
     B1 = B01*diag(w0);
     L0 = B1*B1';
end







