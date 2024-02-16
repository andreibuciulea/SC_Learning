function out = chepuri_alg(X0,X1,prms)
    E = prms.E;
    [N,~] = size(X0);    
    L1 = (ones(N) - diag(ones(1,N))); 
    G1 = graph(L1~=0);
    
    % Incidence matrix of the complete graph
    A1 = full(incidence(G1));
    M =  size(A1,2);
    
    % solution based on sorting
    c=zeros(M,1);
    for ii=1:M
       c(ii) = trace(X'*(A1(:,ii)*A1(:,ii)')*X);
    end
    [~,ind_opt] = sort(c, 'ascend');
    
    w1 = zeros(M,1); w1(ind_opt(1:E)) = 1;
    L_hat = (A1*diag(w1)*A1');
    out.L0 = L_hat;
    H1 = diag(diag(L_hat))-L_hat;
    out.H1 = H1;
    %out.H2 = (generate_B(ones(N)).B0'.*kr(H1,H1))';
    %out.H2kr = out.H2;
    
    prms.L_hat = L_hat; 
    bar_out = barbarossa_alg(X1,prms);
    out.H2 = bar_out.H2;
    out.H2kr = out.H2;
end