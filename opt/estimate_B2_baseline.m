function out = estimate_B2_baseline(X1,w1,prms)
    %B1 --> node-edges incidence matrix
    %B02 --> edges-triangles incidence matrix of a full connected graph
    %a1 --> controls the scale of L0
    %rho--> controls the presence of edges in a triangle

    B02 = prms.B02;
    [E,T] = size(B02);
    a1 = prms.a1;
    rho = prms.rho;
    verbose = prms.verbose;
    C = X1*X1';

    %search for all the triangles in which the observed edges are involved
    w1_o = sum(abs(X1'))~=0; %these are the observed edges
    B2_o = diag(w1_o)*B02;
    w2_p = (sum(abs(B2_o))==3)'; %binary vector with 1 for all possible triangles from observed edges
    
    option = 1;

    if option == 1

        c0 = diag(B02'*C*B02*eye(T));%cost for the smoothness for each triangle
        %auxiliar vector for not considering triangles with a smooth value smaller than 1e-5
        c0 =w2_p.*c0;%use only smoothness associated with all posible triangles from observed signals
        %c1 = w0'*abs(B02)*eye(T);c1 = c1'/max(c1);
        c1 = zeros(size(c0));
        c0 = c0/max(c0);
        c2 = a1*c0+ + 5*not(w2_p); %compute the total cost
        [a,~] = sort(c2);c2_bin = c2<=a(prms.T);c2_bin = c2_bin>=1e-5;
        
        if verbose 
            figure(4)
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
            title('c2 aux')
            subplot(414)
            plot(c2_bin);
            ylim([-0.5,1.5])
            title('c2 bin')
            sgtitle('B2')
        end
        
        if sum(w2_p) <= prms.T %if all possible links are T, then use all of them
            w2 = w2_p;
        else %if not,use the smoothest T triangles
            w2 = c2_bin;
        end
    elseif option == 2
    
        cvx_begin quiet
            variables w1(T,1)
            expressions L1(E,E)
            B2 = B02*diag(w2);
            minimize(a1*norm(B2'*X1,"fro") + ...
                            + rho*(1-w1)'*B02*w2) 
            subject to 
                w2 >= 0;
                w2 <= 1;
                sum(w2) == prms.T;
        cvx_end
        
        [a,~] = sort(w2,"descend");w2 = w2>=a(prms.T);w2 = w2>1e-5;
    end

    B2 = B02*diag(w2);
    L1u = B2*B2';

    out.w2 = w2;
    out.B2 = B2;
    out.L1u = L1u;
end