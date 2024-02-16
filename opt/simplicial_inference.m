function out = simplicial_inference(Y0,Y1,E,T,prms)

    %The number of nodes and edges is known
    [N,M] = size(Y0);
    p1 = M*N;
    %C = Y0*Y0'/p1; 
    %E = size(Y1,1); %this E is always N*(N-1)/2
    
    max_iters = prms.max_iters;
    
    %initialization for B01 and B02 for the fully connected graph
    Bout = gen_B12(N);
    B01 = Bout.B1;B02 = Bout.B2;
    

    %%% B2 related to the fully connected graph
    B1 = B01;B2 = B02;
    T0 = size(B2,2);
    w2 = zeros(T0,1);
   
    %%%these edges should be always there
    

    %%% P0 binary with the entries related to the observed node signals
    P0 = ones(size(Y0));
    p0idx = find(Y0==0); %find the unobserved entries 
    P0(p0idx) = 0;
    X0 = Y0; %Y0 is the original data, X0 is the estimaded node data

    %%% P1 binary with the entries related to the observed edge signals
    P1 = zeros(size(Y1));
    p1idx = find(sum(Y1'~=0)~=0); %prior knowledge about the edges present in the graph
    P1(p1idx,:) = 1;
    X1 = Y1; %Y1 is the original edge data, X1 is the estimaded edge data

    prms.p1idx = p1idx;
    prms.B01 = B01;
    prms.B02 = B02;
    prms.E = E;
    prms.T = T;
    for i = 1:max_iters
        
    %%%%%% Step 1 %%%%%%
        % First iteration we use full B1 to estimate X0
        %X0 = estimate_X0(P0,Y0,B1,prms).X0;
        X0=Y0; %Added by AGM, to be removed
        
    %%%%%% Step 2 %%%%%%
        %%% Estimate B1 from estimated X0 and full B2 for the first iteration 
        C = X0*X0'/p1; 
        B1out = estimate_B1(C,w2,prms);

        w1 = B1out.w0;
        B1 = B1out.B1;
        L0 = B1out.L0;

    %%%%%% Step 3 %%%%%%
        %%% Estimate X1 from full B2 in the first iteration 
        X1 = estimate_X1(P1,Y1,B2,prms).X1;
        %X1=Y1; %Added by AGM, kills the optimization
        
    %%%%%% Step 4 %%%%%%
        %%%Estimate B2 using the estimated B1 and imposing orthogonality
        B2out = estimate_B2(X1,w1,prms);
        w2 = B2out.w1;
        B2 = B2out.B2;
        L1u = B2out.L1u;
        
        %fobj(i) = vec(C)'*vec(L0) + trace(Xsh'*(B02*diag(w1)*B02')*Xsh)/p2;
        %disp(['Error:' num2str(fobj(i)) ' T =' num2str(T)])

%         figure(1)
%         subplot(221)
%         imagesc(prms.B1*prms.B1')
%         colorbar()
%         title('L0 true')
%         subplot(222)
%         imagesc(L0)
%         colorbar()
%         title('L0 hat')
%         subplot(223)
%         imagesc(abs(prms.B1*prms.B1'-L0))
%         colorbar()
%         title('L0true-L0hat')
%         subplot(224)
%         imagesc(L1u)
%         colorbar()
%         title('L1u = B2*B2')
%         sgtitle(['Estimated L0 and L1u at iteration: ' num2str(i)])
    end
    %disp(['Number of estimated triangles: ' num2str(size(B02*diag(w1),2))]) %Original version by AB
    %disp(['Number of estimated triangles: ' num2str(sum(sum(B02*diag(w1))))]) %Added by AGM
    out.B1 = B1;
    out.B2 = B02(:,logical(w2));
    out.X0 = X0; %Added by AGM
    out.X1 = X1; %Added by AGM
    out.w1 = w1;
    out.w2 = w2;
    
end