function out = sim_inf(in)
    Y0 = in.Y0;
    Y1 = in.Y1;
    N = in.N;
    Bout = gen_B12(N);
    B1 = Bout.B1;
    B2 = Bout.B2;
    [E,T] = size(B2);
    max_iters = 10;
    a0 = 1e-3;
    a1 = 1e-3;
    a2 = 1e-3;
    a3 = 1;
    a4 = 1;
    a5 = 1;
    X1 = Y1;
    %Y0 = [Y0;randn(H_0,M_0)]; %complete Y0 with zeros
    P0 = eye(N);%Construct P0 
    Ie = eye(E);
    In = eye(N);
    onET = ones(E,T);
    onNE = ones(N*E,1);
    for i = 1:max_iters
        %%%%%%%%%%%%%optimize w.r.t X0 fixing X1, B1, and B2
        X0 = inv(eye(N)+a0*B1*B1')*Y0;
        
        %%%%%%%%%%%%%optimize w.r.t B1 fixing X0, X1, and B2
        A = a0*(X0*X0');
        C = a1*(X1*X1')+a2*(B2*B2');
        b1 = -inv(kron(Ie,A)+kron(C',In))*a3*onNE;
        b1(b1 < -1) = -1;
        b1(b1 > 1) = 1;
        B1 = reshape(b1,[N,E]);
        %project
        
        %%%%%%%%%%%%%optimize w.r.t X1 fixing X0, B1, and B2
        P1 = diag(ones(E,1));%Construct P1     
        X1 = inv(eye(E)+P1'*P1+a1*B1'*B1+a2*B2*B2')*P1'*P1*Y1;

        %%%%%%%%%%%%%optimize w.r.t B2 fixing X0, X1, and B1
        B2 = -inv(a4*(X1*X1'))*a5*onET;
        B2(B2>1) = 1;
        B2(B2<-1) = -1;

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