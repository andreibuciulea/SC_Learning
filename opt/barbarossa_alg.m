function out = barbarossa_alg(X,prms)
    T = prms.T;%number of triangles
    L_hat = prms.L_hat;
    est_A = diag(diag(L_hat))-L_hat;
    [E,~] = size(X);    
    N = sqrt(E);
    % Incidence matrix of edges vs triangles
    A = ones(N)-eye(N); 
    T = trace(A^3)/6;
    B2 = zeros(E,T);
    R = 1:N;
    t = 1;
    for i = R
        for j = R(R>i)
            for k = R(R>i & R>j)
                %disp(['[' num2str(i) ','  num2str(j) ','  num2str(k) ']'])
                idx(1) = N*(i-1)+j;
                idx(2) = N*(i-1)+k;
                idx(3) = N*(j-1)+k;
                B2(idx(1:2),t) = 1;
                B2(idx(3),t) = -1;
                t = t+1;
            end
        end
    end
    
    %incidence matrix of nodes vs edges
    B1 = zeros(N,E);
    for i = R
        for j = R(R>i)
            if est_A(i,j) ~= 0
                B1(i,N*(i-1)+j) = -1;
                B1(j,N*(i-1)+j) = 1;
            end
        end
    end

    % solution based on sorting
    c=zeros(T,1);
    for ii=1:T
       c(ii) = trace(X1'*(B2(:,ii)*B2(:,ii)')*X1);
    end
    [~,ind_opt] = sort(c, 'ascend');
    
    w1 = zeros(T,1); w1(ind_opt(1:T)) = 1;
    L_hat= (B2*diag(w1)*B2');
    out.LB2 = L_hat;
    H1 = diag(diag(L_hat))-L_hat;
    out.H1 = H1;
    
end