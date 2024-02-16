function out = estimate_B2(X1,w1,prms)
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
    
    %hay que ver que pasa con B
    %mirar los triangulos bien estimados
    %comprobar los triangulos más smooth

    option = 1;

    if option == 1
        %c0=zeros(T,1);
        %for ii=1:T
        %   c0(ii) = trace(X1'*(B02(:,ii)'*B02(:,ii))*X1);
        %end
        c0 = diag(B02'*C*B02*eye(T));%cost for the smoothness for each triangle
        %auxiliar vector for not considering triangles with a smooth value smaller than 1e-5
        c_aux = c0>1e-5;%this is because we do not have access to all the edge signals
        c1 = w1'*abs(B02)*eye(T);
        c0 = c0/max(c0);
        c1 = c1'/max(c1);
        c2 = zeros(size(c0));
        %c2 = a1*c0 - rho*c1;
        c2(c_aux) = a1*c0(c_aux) - rho*c1(c_aux); %compute the total cost
        [~,b] = sort(c2);
        c2_bin = zeros(size(c0));
        c2_bin(b(1:prms.T)) = 1;
        %c2_bin = c2<=a(prms.T);
        %c2_bin = c2_bin>=1e-5;
        
        if verbose 
            figure(8)
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

       w2 = c2_bin;
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

    out.w1 = w2;
    out.B2 = B2;
    out.L1u = L1u;
end