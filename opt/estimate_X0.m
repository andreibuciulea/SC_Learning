function out = estimate_X0(P0,Y0,B1,prms)
    %P0  --> prior information about the missing entries in Y0
    %Y0  --> noisy/incomplete node signal information
    %B1  --> node-edges incidence matrix
    %la0 --> controlling similarity between Y0 and X0 
    %a0  --> controlling node signal smoothness
    la0 = prms.la0;
    a0 = prms.a0;
    [N,M] = size(Y0);

    option = 1;

    if option == 1
        %disp('Close form solution')
        V = la0*diag(P0(:)) + kron(eye(M),a0*B1*B1');%inverse of the sum of matrices
        MatToVec=la0*P0.*Y0;
        x0 = V\MatToVec(:);
        X0 = reshape(x0,[N,M]);
    elseif option == 2
        cvx_begin quiet
            variable X0(N,M)
            minimize(la0*norm(P0.*(Y0-X0),'fro') + a0*norm(B1'*X0,"fro"))
        cvx_end
    end
    out.X0 = X0;
end