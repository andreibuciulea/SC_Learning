function out = est_B1(X0,X1,prms)
    K = prms.K;
    [N,~] = size(X0);    
    L1 = (ones(N) - diag(ones(1,N))); 
    G1 = graph(L1~=0);
    
    % Incidence matrix of the complete graph
    A1 = full(incidence(G1));
    M =  size(A1,2);
    
    % solution based on sorting
    c0=zeros(M,1);
    for ii=1:M
       c0(ii) = trace(X0'*(A1(:,ii)*A1(:,ii)')*X0);
    end
    c1=zeros(N,1);
    for ii=1:N
       c1(ii) = trace(X1'*(A1(ii,:)'*A1(ii,:))*X1);
    end
    [~,ind_opt] = sort(c0, 'ascend');
    
    %es necesario minimizar c1 tambien de alguna manera.
    %se podría seleccionar los enlaces que presentan un 10%
    %del valor de smoothness total cumsum(sort(c,'ascend'))

    w1 = zeros(M,1); w1(ind_opt(1:K)) = 1;
    L_hat = (A1*diag(w1)*A1');
    out.L_hat = L_hat;
    H1 = diag(diag(L_hat))-L_hat;
    out.H1 = H1;
    out.B1 = A1*diag(w1);
end