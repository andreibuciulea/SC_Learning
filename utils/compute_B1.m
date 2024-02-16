function out = compute_B1(A)
    N = size(A,1);
    E = N*(N-1)/2;
    B1 = zeros(N,E);
    NN = 1:N;
    e = 1;
    for a = NN
        for b = NN(NN>a)
            if(A(a,b)~= 0)
                B1(a,e) = -1;
                B1(b,e) = 1;
            end
            e = e+1;
        end
    end
    out.B1 = B1;
    out.w1 = ~all(B1==0)';
end