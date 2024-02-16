function out = lap_inf(X0,X1,prms)
    %The number of nodes and edges is known
    [N,M] = size(X0);
    p1 = M*N;
    C = X0*X0'/p1; 
    E = size(X1,1); %this E is always N*(N-1)/2
    
    a0 = prms.a0;
    a1 = prms.a1;
    a2 = prms.a2;
    
    p2 = M*E;
    max_iters = 30;%prms.max_iters;
    
    %initialization for B01 and B02 for the fully connected graph
    Bout = gen_B12(N);
    B01 = Bout.B1;B02 = Bout.B2;

    idx = []; %no prior knowledge about the edges present in the graph
    for i = 1:max_iters
        cvx_begin quiet
            variables w0(E,1)
            expressions L0(N,N) symmetric semidefinite
            L0 = B01*diag(w0)*B01';
            minimize(a0*vec(C)'*vec(L0) + a1*vec(L0(~eye(N)))'*vec(L0(~eye(N)))...
                - ones(1,N)*log(diag(L0)))
            subject to
                %sum(abs(L0*ones(N,1))) <= 1e-8;
                %L0(~eye(N)) <= 0;
                w0 >= 0;
                w0(idx) == 1;
         cvx_end

        w0 = w0>1e-5;
        %update B1 
        B1 = B01*diag(w0);
        L0 = B1*B1';
        
        A = diag(diag(L0))-L0;
        %update B02 (with the new value for T)
        B02 = compute_B2(A,B1);
        %update T
        T = size(B02,2);
        %project the signal
        [V1,V2] = eig(B1'*B1);
        cols = find(diag(V2)> 1e-5);
        Uirr = V1(:,cols);
        Xsh = (eye(E)-Uirr*Uirr')*X1;

        cvx_begin quiet
            variables w1(T,1)
            expressions L1(E,E)
            L1 = B02*diag(w1)*B02';
            minimize(trace(Xsh'*(B02*diag(w1)*B02')*Xsh)/p2 + a2*norm(L1,"fro"))
         
            subject to 
                w1 >= 0;
                sum(w1) == 1;
                %trace(L1) >= N; %L1 == L1';%L0(~I) <= 0;%trace(L0) >= N;
        cvx_end
        
        %binarize w1?
        %w1 = w1>prctile(w1,5);
        w1 = w1>1e-5;
        %figure(1)
        %plot(w1)
        %update idx with the links present in triangles
        Aux = B02*diag(w1);
        idx = sum(Aux,2) ~= 0;

        fobj(i) = vec(C)'*vec(L0) + trace(Xsh'*(B02*diag(w1)*B02')*Xsh)/p2;
        disp(['Error:' num2str(fobj(i)) ' T =' num2str(T)])

%         figure(1)
%         subplot(131)
%         imagesc(L0)
%         colorbar()
%         subplot(132)
%         imagesc(L0)
%         colorbar()
%         subplot(133)
%         imagesc(L1)
%         colorbar()
%         title(' ')
    end
    disp(['Number of estimated triangles: ' num2str(size(B02*diag(w1),2))])
    out.B1 = B1;
    out.B2 = B02*diag(w1);

end