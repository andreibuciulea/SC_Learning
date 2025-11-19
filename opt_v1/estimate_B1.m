function out = estimate_B1(X0,w2,prms)
%Input
    %X0 --> noisy node data
    %w1 --> triangle selection vector
    %B01 --> node-edges incidence matrix of fully connected graph
    %B02 --> edges-triangles incidence matrix of a fully connected graph
    %prms --> structure containing parameters
        %alpha1 --> controls the norm of w1
        %beta1 --> controls the smoohness of X0
        %gamma --> controls valid SC
        %E --> number of edges
        %p1idx --> prior knowledge about the existing edges
%Output

    N = size(X0,1);

    alpha1 = prms.alpha1;
    beta1 = prms.beta1;
    gamma = prms.gamma;
    p1idx = prms.p1idx;
    verbose = prms.verbose;
    B01 = prms.B01;
    B02 = prms.B02;
    B01t = B01';
    w1_true = prms.w1;% this is the true w1
    
    option = 1;
    %some normalization factor should be added here:X0*X0'
    if option == 1
        t0 = tic();
        %c0 = diag(B01'*(X0*X0')*B01);% cost for the smoothness for each edge
        Y = B01t * X0; 
        c0 = sum(Y.^2, 2);
        c0 = c0/max(max(c0));
        c1 = abs(B02)*w2;%cost for the presence of an edge in the triangle
        c1 = c1/max(max(c1));
        s1 = alpha1 + beta1*c0 - gamma*c1; %compute the total cost
        %The values associated with prior information should be -1 
        s1(p1idx) = -1;
        %s1 can be positive and negative 
        %We should select all edges with negative cost in s1 (no matter if the number of edges is larger than prms.E)
            %if the number of these entries is less than prms.E
            %then we select smallest positive entries until we reach prms.E
        % Sort the elements of s1 in ascending order and get their indices
        [~, b] = sort(s1);
        
        % Initialize the binary selection vector with zeros
        s2_bin = zeros(size(c0));
        % Select all edges with negative cost
        negative_indices = find(s1 < 0);
        % Determine the number of edges to be selected
        num_negative = length(negative_indices);
        E = prms.E;  % Number of edges to be selected

        if num_negative >= E
            % If there are at least E negative edges, select the E smallest ones
            s2_bin(b(1:E)) = 1;
        else
            % If there are fewer than E negative edges, select all of them
            s2_bin(negative_indices) = 1;
            % Select additional smallest positive entries to reach E edges
            remaining = E - num_negative;
            positive_indices = find(s1 >= 0);
            % Sort positive entries to get the smallest ones
            [~, pos_sort_idx] = sort(s1(positive_indices));
            % Select the smallest positive edges needed
            s2_bin(positive_indices(pos_sort_idx(1:remaining))) = 1;
        end
        w1 = s2_bin;
        t1 = toc(t0);
    elseif option == 2
        cvx_begin quiet
            variables w0(E,1)
            expressions L0(N,N) symmetric semidefinite
            L0 = B01*diag(w1)*B01';
            minimize(a0*vec(C)'*vec(L0) + ga*vec(L0(~eye(N)))'*vec(L0(~eye(N)))...
                - ones(1,N)*log(diag(L0)) + rho*(1-w1)'*B02*w1)
            subject to
                w1 >= 0;
                w1(p1idx) == 1; %prior edge knowledge
        cvx_end
    
        %%% Update B1
        [a,~] = sort(w1,"descend");w1 = w1>=a(prms.E);w1 = w1>=1e-5;
    end
    
    [~,r_e] = ranking_edges(b,w1,w1_true);

    B1 = B01*diag(w1);
    L0 = B1*B1';
    
    out.w1 = w1;
    out.B1 = B1;
    out.L0 = L0;
    out.re = r_e;
    out.t1 = t1;
end