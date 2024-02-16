function out = estimate_SC(X0,X1,E,T,model,prms,B2,w1)
    if strcmp(model,'ours')
        out = simplicial_inference(X0,X1,E,T,prms);
    elseif strcmp(model,'baseline')
        out = baseline_SC(X0,X1,E,T,prms);
    elseif strcmp(model,'check')
        out = check(X0,X1,B2,w1);
    elseif strcmp(model,'barbarossa')
        disp('Not available yet')
    elseif strcmp(model,'RC')
        disp('Not available yet')
    else
        disp('Unknown estimation method')
    end
end