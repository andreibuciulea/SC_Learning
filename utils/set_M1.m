function wm1 = set_M1(w1,ps)

    E = sum(w1);
    nl = round(ps*E);
    [a,~] = find(w1==1);
    idx = randperm(E,nl);
    wm1 = zeros(size(w1));
    wm1(a(idx)) = 1; %these are the selected edges

end