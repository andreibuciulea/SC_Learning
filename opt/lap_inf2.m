function out = lap_inf2(Y0,Y1,prms)
    %hay que penalizar de manera distinta un link que se añade para formar
    %un triangulo que un link que se quita para romper un triangulo.
    %Y también diferente a un link que se añade o quita sin tener que ver
    %con los triangulos.
    
    %X0 y X1 no están completos.

    %The number of nodes and edges is known
    [N,M] = size(Y0);
    p1 = M*N;
    C = Y0*Y0'/p1; 
    E = size(Y1,1); %this E is always N*(N-1)/2
    
    a0 = prms.a0;
    a1 = prms.a1;
    a2 = prms.a2;
    rho = 1e-3;
    
    p2 = M*E;
    max_iters = 2;%prms.max_iters;
    
    %initialization for B01 and B02 for the fully connected graph
    Bout = gen_B12(N);
    B01 = Bout.B1;B02 = Bout.B2;

    %%% B2 related to the fully connected graph
    B1 = B01;B2 = B02;
    T = size(B2,2);
   
    %%%these edges should be always there
    

    %%% P0 binary with the entries related to the observed node signals
    P0 = ones(size(Y0));
    p0idx = find(Y0==0); %find the unobserved entries 
    P0(p0idx) = 0;
    X0 = Y0; %Y0 is the original data, X0 is the estimaded node data

    %%% P1 binary with the entries related to the observed edge signals
    P1 = ones(size(Y1));
    p1idx = [];%find(sum(Y1'~=0)~=0); %prior knowledge about the edges present in the graph
    P1(p1idx) = 0;
    X1 = Y1; %Y1 is the original edge data, X1 is the estimaded edge data

    for i = 1:max_iters

        %%%%%% Step 1 %%%%%%
        %%% first iteration we use full B1 to estimate X0
        la0 = 1; %this parameter should be used to control similarity between X0 and Y0
        cvx_begin quiet
            variable X0(N,M)
            minimize(la0*norm(P0.*(Y0-X0),'fro') + a0*norm(B1'*X0,"fro"))
        cvx_end

        X0=Y0; %Added by AGM, to me removed
        %%%%%% Step 2 %%%%%%
        %%% Estimate B1 from estimated X0 and full B2 for the first iteration 
        C = X0*X0'/p1; 
        cvx_begin quiet
            variables w0(E,1)
            expressions L0(N,N) symmetric semidefinite
            L0 = B01*diag(w0)*B01';
            minimize(a0*vec(C)'*vec(L0) + a1*vec(L0(~eye(N)))'*vec(L0(~eye(N)))...
                - ones(1,N)*log(diag(L0)) + rho*norm((B01*diag(w0))*B2,"fro"))
            subject to
                %sum(abs(L0*ones(N,1))) <= 1e-8;
                %L0(~eye(N)) <= 0;
                w0 >= 0;
                w0(p1idx) == 1; %prior edge knowledge
        cvx_end
        %%% Update B1
        [a,~] = sort(w0,"descend");w0 = w0>=a(prms.E);w0 = w0>=1e-5;
        
        B1 = B01*diag(w0);
        B1=prms.B1; %Added by AGM, to me removed

        L0 = B1*B1';
        
        
        
        %%%%%% Step 3 %%%%%%
        %%% Estimate X1 from full B2 in the first iteration 
        %%% Proyect signal X1
        [V1,V2] = eig(B1'*B1);cols = find(diag(V2)> 1e-5);
        Uirr = V1(:,cols);
        H = (eye(E)-Uirr*Uirr');
        la1 = 1; %this parameter should be used to control similarity between X1 and Y1
        cvx_begin quiet
            variable X1(E,M)
            minimize(la1*norm(P1.*(Y1-X1),'fro') + a1*norm(B2'*(H*X1),"fro"))
        cvx_end
        
        X1=Y1; %Added by AGM, kills the optimization
        
        %%%%%% Step 4 %%%%%%
        %%%Estimate B2 using the estimated B1 and imposing orthogonality
        Xsh = H*X1;
        Xsh = X1; %Added by AGM kills the previous line
        penalty_triang = compute_penalty_triangles(B1,B02); %Created by AGM
        cvx_begin quiet
            variables w1(T,1)
            expressions L1(E,E)
            L1 = B02*diag(w1)*B02';
%            minimize(1e5*trace(Xsh'*(B02*diag(w1)*B02')*Xsh)/p2 + a2*norm(L1,"fro")...
%                            + rho*norm(B1*(B02*diag(w1)),"fro"))
            minimize(1e0*trace(Xsh'*(B02*diag(w1)*B02')*Xsh)/p2 + ...
                            + rho*trace(B1*(B02*diag(w1)*B02'*B1')) + ... 
                        + rho*penalty_triang*w1)%AGM: replaced the objective 

                        %Here we have a problem because
                        %diag(B02'*B1'*B1*B02) is zero for triangles
                        %involving nonexisting links in B1 ==> for this
                        %simple (sorted formulation), we can use function
                        %compute penalty_triangles)
            subject to 
                w1 >= 0;
                % sum(w1) == 1; %AGM commented this
                w1 <= 1;
                sum(w1) == prms.T;
                %trace(L1) >= N; %L1 == L1';%L0(~I) <= 0;%trace(L0) >= N;
        cvx_end
        
        [a,~] = sort(w1,"descend");w1 = w1>=a(prms.T);%w1 = w1>1e-5;
        B2 = B02*diag(w1);
        
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
    disp(['Number of estimated triangles: ' num2str(size(B02*diag(w1),2))]) %Original version by AB
    disp(['Number of estimated triangles: ' num2str(sum(sum(B02*diag(w1))))]) %Added by AGM
    out.B1 = B1;
    out.B2 = B02*diag(w1);
    out.X0 = X0; %Added by AGM
    out.X1 = X1; %Added by AGM
    

end