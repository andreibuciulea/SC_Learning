function out = check(X0,Xm1,B2,w1_true,w2_true)
    
    %Is obtaining the smoothest graph without considering knowledge about
    %the observed edges, only node signals
    %Is obtaining the smoothest triangles without considering knowledge
    %about the graph connectivity, only X1.
    
    N = size(X0,1);
    Bout = gen_B12(N);
    B01 = Bout.B1;B02 = Bout.B2;
    Theta = diag(w1_true);

    [w2_s,b2] = check_smooth_X1(Xm1,B2,N);
    [w1_s,b1] = check_smooth_X0(X0,w1_true,N);

    [~,r_t] = ranking_triangles(b2,w2_s,w2_true);
    [~,r_e] = ranking_edges(b1,w1_s,w1_true);

    out.w1 = w1_s;
    out.w2 = w2_s;
    B1 = B01*diag(w1_s);
    out.B1 = B1;
    B2 = B02(:,w2_s);
    out.B2 = B2;
    out.X1 = Xm1;
    out.X0 = X0;
    out.rt = r_t;
    out.re = r_e;
    %out.X1 = pinv(Theta'*Theta + (B2*B2'))*(Theta'*Xm1);
    %out.X0 = (eye(N) + (B1*B1'))\X0;
end