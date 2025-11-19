function out = estimate_B2(X1,w1,prms)
%Input
    %X1 --> noisy observed edge data
    %w1 --> edge selection vector
    %B02 --> edges-triangles incidence matrix of a fully connected graph
    %prms --> structure containing parameters
        %alpha2 --> controls the norm of w2
        %beta2 --> controls the smoohness of X1
        %gamma --> controls valid SC
        %T --> number of triangles
%Output
    M = size(X1,2);
    B01 = prms.B01;
    B02 = prms.B02;
    T = size(B02,2);
    alpha2 = prms.alpha2;
    beta2 = prms.beta2;
    gamma = prms.gamma;
    verbose = prms.verbose;
    w2_true = prms.w2;
    B02t = B02';

    option = 1;
    %some normalization factor should be added here:X1*X1'
    if option == 1
        t0 = tic();
        %c0 = diag(B02'*(X1*X1')*B02)/M;%cost for the smoothness for each triangle
        Y = B02t * X1;       
        c0 = sum(Y.^2, 2);  
        c0 = c0/max(max(c0));
        c1 = abs(B02)'*(1-w1); %cost for a valid SC 
        s2 = alpha2 + beta2*c0 + gamma*c1; %compute the total cost
        %b are the indexes of the estimated triangles sorted
        [a,b] = sort(s2);
        s2_bin = zeros(size(c0));
        s2_bin(b(1:prms.T)) = 1;
        w2 = s2_bin;
        t1 = toc(t0);
    elseif option == 2
        % This should be revised
        cvx_begin quiet
            variables w1(T,1)
            expressions L1(E,E)
            B2 = B02*diag(w2);
            minimize(a1*norm(B2'*X1,"fro") + ...
                            + rho*(1-w2)'*B02*w2) 
            subject to 
                w2 >= 0;
                w2 <= 1;
                sum(w2) == prms.T;
        cvx_end
        
        [a,~] = sort(w2,"descend");w2 = w2>=a(prms.T);w2 = w2>1e-5;
    end


    [~,r_t] = ranking_triangles(b,s2_bin,w2_true);

    B2 = B02*diag(w2);
    L1u = B2*B2';
    w1T = sum(abs(B2),2)~=0;
    w1(w1T) = 1;
    B1 = B01*diag(w1);
    out.w1 = w1; %these are the edges that should be present for the filled triangles to exist 
    out.w2 = w2;
    out.B1 = B1;
    out.B2 = B2;
    out.rt = r_t;% average percentage of the ranking for the incorrectly estimated triangles
                %if you have 10% and 1000 possibles triangles, the
                %incorrectly estimated triangles, on average, are in the
                %first 100 positions after the true ones.

    out.L1u = L1u;
    out.t1 = t1;
end
