function out = estimate_X0(P0, Y0, B1, prms)
    %P0  --> prior information about the missing entries in Y0
    %Y0  --> noisy/incomplete node signal information
    %B1  --> node-edges incidence matrix
    %prms --> structure containing parameters
        % --> beta1: smoothness regularizer for X0 
        % --> eta0: noise regularizer for X0
        % --> max_iters: for iterative approaches
        % --> verbose: true for print results
    
    % Default values for parameters
    %this should be done properly by checking if the parameter is there if
    %not default value is assigned
    default_values = struct('beta1', 1, 'eta0', 1, 'max_iters', 100,'verbose',false);
    
    % Check and assign parameters
    if prms.verbose
        param_names = fieldnames(default_values);
        for i = 1:numel(param_names)
            if ~isfield(prms, param_names{i})
                prms.(param_names{i}) = default_values.(param_names{i});
                disp(['X0 estimation: Parameter ', param_names{i}, ' is not provided. Using default value: ', num2str(default_values.(param_names{i}))]);
            end
        end
    end
    beta1 = prms.beta1;
    eta0 = prms.eta0;
    max_iters = prms.max_iters;
    verbose = prms.verbose;
    
    [N, M] = size(Y0);

    %this should be a parameter
    option = 1;

    %This should be revised: optimal choice for mu
    %lambdas = eig(2*a0*(B1*B1') + 2*la0*eye(N));
    %mu = 1/max(lambdas);
    if option == 1
        X0 = (eye(N) + beta1/eta0*(B1*B1'))\Y0; 
    elseif option == 2
        % This should be revised
        V = la0 * diag(P0(:)) + kron(eye(M), a0*(B1*B1')); % inverse of the sum of matrices
        MatToVec = la0 * P0 .* Y0;
        x0 = V \ MatToVec(:);
        X0 = reshape(x0, [N, M]);
    elseif option == 3
        % This should be revised
        cvx_begin quiet
            variable X0(N, M)
            minimize(la0 * norm(P0 .* (Y0 - X0), 'fro') + a0 * norm(B1' * X0, "fro"))
        cvx_end
    elseif option == 4
        % This should be revised
        X = Y0;
        for i = 1:max_iters
            X0 = X-2*mu*(a0*(B1*B1')*X+la0*(Y0-X));
            xn = norm(X-X0,'fro');
            X = X0;
            if xn <= th
                if verbose
                    disp(['Converged in ' num2str(i) ' iterations']);
                end
                break;
            end
        end
    end
    
    out.X0 = X0;
end
