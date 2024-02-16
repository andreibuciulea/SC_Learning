function out = select_hidden_nodes(L,X,H,sel_type)
    N = size(L,1);
    if nargin < 3
       H = 1;
       sel_type = 'random';
    end

    switch sel_type
        case 'random'
            v = randperm(N,N);
        case 'min_d'
            [~,v] = sort(diag(L));
        case 'max_d'
            [~,v] = sort(diag(L),'descend');
        otherwise
            error('ERR: Unkown selection type')
    end
    s_h = sort(v(1:H));
    s_o = sort(v(H+1:end));
    
    out.Lo = L(s_o,s_o);
    out.Xo = X(s_o,:);
    out.s_o = s_o;
    out.s_h = s_h;
end