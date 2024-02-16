function out = est_B2(in)
    X1 = in.X1; 
    K = in.T;%number of triangles
    L = in.L;
    A = diag(diag(L))-L;
    [E,~] = size(X1);    
    
    %incidence matrix of nodes vs edges
    B1 = compute_B1(A).B1;

    % Incidence matrix of edges vs triangles
    B2 = compute_B2(A,B1);
    T = size(B2,2);
    %use B1 to remove the irrotational part from X
    [U,Sig] = eig(B1'*B1);
    %select the eigenvector associated to nonzero eigenvalues 
    idx = find(diag(Sig)>1e-6);
    Uirr = U(:,idx);
    Xsh = (eye(E)-Uirr*Uirr')*X1;
    % solution based on sorting
    c=zeros(T,1);
    for ii=1:T
       c(ii) = trace(Xsh'*(B2(:,ii)*B2(:,ii)')*Xsh);
    end
    [~,ind_opt] = sort(c, 'ascend');
    
    w1 = zeros(T,1); w1(ind_opt(1:K)) = 1;
    L_hat= (B2*diag(w1)*B2');
    out.L_hat = L_hat;
    H1 = diag(diag(L_hat))-L_hat;

    out.H1 = H1;
    out.B2 = B2*diag(w1);
end