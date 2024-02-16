function out = simplicial_inf(in)
    Y0 = in.Y0;
    Y1 = in.Y1;
    N = in.N;
    w1 = in.w1;
    [O_0,M_0] = size(Y0);
    [O_1,M_1] = size(Y1);
    H_0 = N-O_0;
    E = N*(N-1)/2;
    Bout = gen_B12(N);
    B_01 = Bout.B1;
    B_02 = Bout.B2;
    B1 = B_01;
    B2 = B_02;
    max_iters = 10;
    a0 = 1e-2;
    a1 = 10;
    a2 = 10;
    a3 = 1e-1;
    X1 = Y1;
    Y0 = [Y0;randn(H_0,M_0)]; %complete Y0 with zeros
    P0 = zeros(N);P0(1:O_0,1:O_0) = eye(O_0);%Construct P0 
    for i = 1:max_iters
        %%%%%%%%%%%%%optimize w.r.t X0 fixing X1, B1, and B2
        %start with B1(NxE) as fully connected graph
        %cvx_begin quiet
        %    variable X0(N,M_0)
        %    minimize(norm(P0*(Y0-X0),'fro') + a0*norm(B1'*X0,'fro'))
        %cvx_end
        X0 = inv(eye(N)+a0*B1*B1')*Y0;
        
        %%%%%%%%%%%%%optimize w.r.t B1 fixing X0, X1, and B2
%         cvx_begin quiet
%             variable w(E,1)
%             variable B1(N,E)
%             B1 == B_01*diag(w);
%             minimize(a0*norm(B1'*X0,'fro') + a1*norm(B1*X1,'fro')+ 1e3*norm(w,1))
%             subject to 
%                 B1*B2 == 0;
%                 w >= 0;
%                 sum(w) >= N;
%         cvx_end
        
        K = 39; %do not fix this value
        params = struct('K',K);
        B1 = est_B1(X0,Y1,params).B1;

      
        %%%%%%%%%%%%%optimize w.r.t X1 fixing X0, B1, and B2
        P1 = diag(w1);%Construct P1 
        %cvx_begin quiet
        %    variable X1(E,M_1)
        %    minimize(norm(P1*(Y1-X1),'fro') + a1*norm(B1*X1,'fro') + a2*norm(B2'*X1,'fro') + a3*norm(X1,'fro'))
        %cvx_end      
        X1 = inv(eye(E)+P1'*P1+a1*B1'*B1+a2*B2*B2')*P1'*P1*Y1;
        %%%%%%%%%%%%%optimize w.r.t B2 fixing X0, X1, and B1
        T = 6;
        input = struct('X1',X1,'L',B1*B1','T',T);
        B2 = est_B2(input).B2;

        figure(3)
        subplot(221)
        imagesc(B1)
        colorbar()
        subplot(222)
        imagesc(X0)
        colorbar()
        subplot(223)
        imagesc(B2)
        colorbar()
        subplot(224)
        imagesc(X1)
        colorbar()
    
    end
    out.X0 = X0;
    out.X1 = X0;
    out.B1 = B1;
    out.B2 = X0;
end