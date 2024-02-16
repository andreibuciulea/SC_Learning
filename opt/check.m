function out = check(X0,Xm1,B2,w1)
    
    %Is obtaining the smoothest graph without considering knowledge about
    %the observed edges, only node signals
    %Is obtaining the smoothest triangles without considering knowledge
    %about the graph connectivity, only X1.
    
    N = size(X0,1);
    Bout = gen_B12(N);
    B01 = Bout.B1;B02 = Bout.B2;

    w2_s = check_smooth_X1(Xm1,B2,N);
    w1_s = check_smooth_X0(X0,w1,N);

    out.w1 = w1_s;
    out.w2 = w2_s;
    out.B1 = B01*diag(w1_s);
    out.B2 = B02(:,w2_s);
    out.X1 = Xm1;


end