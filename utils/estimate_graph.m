function out = estimate_graph(X,alg_type,prms)
    
    switch alg_type
        case 'gti_alg'
            gti = gti_alg(X,prms);
        otherwise
            error('ERR: Unkown graph estimation type')
    end
    out.L_hat = gti.L_hat;
end