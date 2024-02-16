function out = estimate_B1(C,w1,prms)
    %C --> sample covariance amtrix of the data
    %B01 --> node-edges incidence matrix of fully connected graph
    %B2 --> edges-triangles incidence matrix (will be removed)
    %a0 --> controlls node signal smoothness
    %a1 --> controls the scale of L0
    %rho--> controls the presence of edges in a triangle
    %p1idx --> prior knowledge about the existing edges

    N = size(C,1);
    a0 = prms.a0;
    ga = prms.ga;
    rho = prms.rho;
    p1idx = prms.p1idx;
    verbose = prms.verbose;
    B01 = prms.B01;
    B02 = abs(prms.B02);
    E = size(B02,1);
    
    option = 1;

    if option == 1
        c0 = diag(B01'*C*B01*eye(E));%cost for the smoothness for each edge
        c1 = eye(E)*B02*w1;%cost for the presence of an edge in the triangle
        c0 = c0/max(c0);
        c1 = c1/max(c1);
        if sum(isnan(c1)) > 0
            c1 = 0;
        end
        %prior information
        cp = zeros(size(c0));
        cp(p1idx) = 1;
        
        c2 = a0*c0 - rho*c1-cp; %compute the total cost
        %select only E edges, no more than that
        [a,b] = sort(c2);
        c2_bin = zeros(size(c0));
        c2_bin(b(1:prms.E))=1;
        %c2_bin = c2<=a(prms.E);
        %c2_bin = c2_bin>=1e-5;%select the E edges with the minimum cost
        w0 = c2_bin;
        if verbose 
            figure(3)
            subplot(411)
            plot(c0);
            ylim([-0.5,1.5])
            title('c0')
            subplot(412)
            plot(c1);
            ylim([-0.5,1.5])
            title('c1')
            subplot(413)
            plot(c2);
            title('c2')
            subplot(414)
            plot(c2_bin);
            ylim([-0.5,1.5])
            title('c2 bin')
            sgtitle('B1')
        end

    elseif option == 2
        cvx_begin quiet
            variables w0(E,1)
            expressions L0(N,N) symmetric semidefinite
            L0 = B01*diag(w0)*B01';
            minimize(a0*vec(C)'*vec(L0) + ga*vec(L0(~eye(N)))'*vec(L0(~eye(N)))...
                - ones(1,N)*log(diag(L0)) + rho*(1-w0)'*B02*w1)
            subject to
                w0 >= 0;
                w0(p1idx) == 1; %prior edge knowledge
        cvx_end
    
        %%% Update B1
        [a,~] = sort(w0,"descend");w0 = w0>=a(prms.E);w0 = w0>=1e-5;
    end
    
    B1 = B01*diag(w0);
    L0 = B1*B1';
    
    out.w0 = w0;
    out.B1 = B1;
    out.L0 = L0;
end