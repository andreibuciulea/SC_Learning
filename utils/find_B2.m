function B2 = find_B2(T,Tv)
    N = 20;
    outB12 = gen_B12(N);
    B01 = outB12.B1;
    B02 = outB12.B2;
    B2 = zeros(size(B02));
    for t = 1:1350
        triangle = T(:,t);
        if sum(triangle) == 0
            break
        end
        %find the links associated with the triangle
        links_idx = find(sum(abs(B01(triangle,:)))==2);
        %buscar el triangulo con los links
        triang_idx = find(sum(abs(B02(links_idx,:)))==3);
        B2(:,triang_idx) = B02(:,triang_idx)*Tv(t);

    end

end