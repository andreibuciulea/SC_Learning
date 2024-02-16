function out = estimate_X1(P1,Y1,B2,prms)
    %P1  --> prior information about the missing entries in Y1
    %Y1  --> noisy/incomplete edge signal information
    %H1  --> proyection matrix for X1
    %B2  --> edge-triangles incidence matrix
    %la1 --> controlling similarity between Y1 and X1 
    %a1  --> controlling edge signal smoothness (cyclic flows)
    [E,M] = size(Y1);
    a1 = prms.a1;
    la1 = prms.la1;

    option = 2;

    if option == 1
        %disp('Close form solution')
        B = la1*diag(P1(:));
        A = kron(eye(M),a1*B2*B2');
        V = A + B;%inverse of the sum of matrices
        MatToVec=la1*P1.*Y1;
        
        V = V + 1e-3*eye(E*M);
        x1 = V\MatToVec(:);

        %V = pinv(V);
        %x1 = V*MatToVec(:);
        
        X1 = reshape(x1,[E,M]);
        %we only have edge signals associated with existing(estimated) edges
        %w1 = sum(abs(B2),2)~=0;
        %X1 = diag(w1)*X1;
    elseif option == 2
        cvx_begin quiet
            variable X1(E,M)
            minimize(la1*square_pos(norm(P1.*(Y1-X1),'fro')) + a1*square_pos(norm(B2'*X1,"fro")))
        cvx_end
    end
    out.X1 = X1;
end