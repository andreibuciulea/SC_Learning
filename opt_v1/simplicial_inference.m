function out = simplicial_inference(Y0,Y1,E,T,prms,B1t,B2t,X0_true,X1t,w1_true,w2_true)
    %Y0 --> noisy node signals (maybe incomplete)
    %Y1 --> noisy edge signals (maybe incomplete, only some of the edges have signals associated)
    %E --> True value for the number of edges
    %T --> True value for the number of triangles
    %prms --> number of iterations

    %The number of nodes and edges is known
    [N,M] = size(Y0);
    p1 = M*N;
    %C = Y0*Y0'/p1; 
    %E = size(Y1,1); %this E is always N*(N-1)/2
    
    max_iters = prms.max_iters;
    %B1t = prms.B1;
    %B2t = prms.B2;
    %X1t = prms.X1;
    
    %initialization for B01 and B02 for the fully connected graph
    Bout = gen_B12(N);
    B01 = Bout.B1;B02 = Bout.B2;
    
    %%% B2 related to the fully connected graph
    B1 = B01;B2 = B02;
    T0 = size(B2,2);
    w2 = zeros(T0,1);
   
    %%%these edges should be always there
    L_true = B1t*B1t';
    S_true = diag(diag(L_true))-L_true;

    %%% P0 binary with the entries related to the observed node signals
    P0 = ones(size(Y0));
    p0idx = find(Y0==0); %find the unobserved entries 
    P0(p0idx) = 0;
    X0 = Y0; %Y0 is the original data, X0 is the estimaded node data
    X1 = Y1;

    %%% P1 binary with the entries related to the observed edge signals
    p1idx = find(sum(abs(Y1')~=0)~=0); %prior knowledge about the edges present in the graph
    % Construct the selection matrix P
    E0 = N*(N-1)/2;
    w1 = zeros(E0,1);
    w1(p1idx) = 1;
    P1 = diag(w1);
    %P1 = eye(size(Y1,1));  % Start with identity matrix
    %P1 = P1(p1idx, :);  % Select the rows corresponding to 'w'
    %Ym1 = Y1(p1idx,:);

    prms.p1idx = p1idx;
    prms.B01 = B01;
    prms.B02 = B02;
    prms.E = E;
    prms.T = T;
    fobj  = zeros(max_iters,3);
    erri  = zeros(4,max_iters);
    all_times = zeros(4,max_iters);
    S0_prev = zeros(N);
    S1_prev = zeros(E0); 
    for i = 1:max_iters
        
    %%%%%% Step 4 %%%%%%
        %%%Estimate B2 using the estimated B1 and imposing orthogonality
        prms.w2 = w2_true;
        B2out = estimate_B2(X1,w1,prms);
        all_times(1,i) = B2out.t1;
        w1 = B2out.w1;
        w2 = B2out.w2;
        B1 = B2out.B1;
        B2 = B2out.B2;
        %L1u = B2out.L1u;

    %%%%%% Step 2 %%%%%%
    %%% Estimate B1 from estimated X0 and full B2 for the first iteration 
        prms.w1 = w1_true;
        t2 = tic();
        B1out = estimate_B1(X0,w2,prms);
        all_times(2,i) = B1out.t1;
        w1 = B1out.w1;
        B1 = B1out.B1;

    %%%%%% Step 1 %%%%%%
        % First iteration we use full B1 to estimate X0
        t3 = tic();
        X0 = estimate_X0(P0,Y0,B1,prms).X0;
        all_times(3,i) = toc(t3);
        %X0=Y0; %Added by AGM, to be removed
        
    %%%%%% Step 3 %%%%%%
        %%% Estimate X1 from full B2 in the first iteration
        t4 = tic();
        X1 = estimate_X1(P1,Y1,B2,prms).X1;
        all_times(4,i) = toc(t4);
        %X1=Y1; %Added by AGM, kills the optimization
        
        fobj(i,1) = compute_obj_function(w1,w2,B01,B02,X0,X1,X0_true,X1t,prms);
      
        L_hat = B1*B1';
        S_hat = diag(diag(L_hat))-L_hat; 

        L1_hat = B2*B2';
        S1_hat = diag(diag(L1_hat))-L1_hat;



        erri(1,i) = norm(B1t-B1,'fro')^2/norm(B1t,'fro')^2;
        erri(2,i) = norm(B2t-B2,'fro')^2/norm(B2t,'fro')^2;
        erri(3,i) = norm(X1t-X1,'fro')^2/norm(X1t,'fro')^2;
        erri(4,i) = norm(S_hat-S_true,'fro')^2/norm(S_true,'fro')^2;
        erri(5,i) = fscore(S_hat,S_true);

        fobj(i,2) = norm(S_hat-S0_prev,'fro')^2;
        fobj(i,3) = norm(S1_hat-S1_prev,'fro')^2;
        S0_prev = S_hat;
        S1_prev = S1_hat;
        

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
    out.obj = fobj;
    out.err = erri;
    out.rt = B2out.rt;
    out.re = B1out.re;
    out.time = all_times;
end