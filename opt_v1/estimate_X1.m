function out = estimate_X1(Theta,Y1, B2, prms)
    %Theta  --> prior information about the missing entries in Y1
    %Y1  --> noisy/incomplete edge signal information
    %B2  --> edge-triangles incidence matrix
    %prms --> structure containing parameters
        % --> beta2: smoothness regularizer for X1 
        % --> eta1: noise regularizer for X1
        % --> max_iters: for iterative approaches
        % --> verbose: true for print results
    
    % Default values for parameters
    default_values = struct('beta2', 1, 'eta1', 1, 'max_iters', 100,'verbose',false);
    
    % Check and assign parameters
    if prms.verbose
        param_names = fieldnames(default_values);
        for i = 1:numel(param_names)
            if ~isfield(prms, param_names{i})
                prms.(param_names{i}) = default_values.(param_names{i});
                disp(['X1 estimation: Parameter ', param_names{i}, ' is not provided. Using default value: ', num2str(default_values.(param_names{i}))]);
            end
        end
    end
    
    beta2 = prms.beta2;
    eta1 = prms.eta1;
    mu1 = 1e-3;%gradient step
    th = 1e-2;
    max_iters = prms.max_iters;
    verbose = prms.verbose;
    
    [E, M] = size(Y1);
    %this should be a parameter
    option = 3;

    %This should be revised: optimal choice for mu
    %lambdas = eig(2*a1*(B2*B2') + 2*la1*(P2'*P2));
    %mu = 1/max(lambdas);
    
    if option == 1
        X1 = pinv(Theta'*Theta + beta2/eta1*(B2*B2'))*(Theta'*Y1);

    elseif option == 4
        % This should be revised
        % Close-form solution
        B = la1 * diag(P1(:));
        A = kron(eye(M), a1 * B2 * B2');
        V = A + B; % inverse of the sum of matrices
        MatToVec = la1 * P1 .* Y1;
        
        V = V + 1e-3 * eye(E * M);
        x1 = V \ MatToVec(:);
        X1 = reshape(x1, [E, M]);
        
    elseif option == 2
        % This should be revised
        cvx_begin quiet
            variable X1(E, M)
            minimize(la1 * square_pos(norm(P1 .* (Y1 - X1), 'fro')) + a1 * square_pos(norm(B2' * X1, "fro")))
        cvx_end
    
    elseif option == 3
        % This should be revised
        X = Y1;
        A = (eye(E)-mu1*beta2*(B2*B2')-mu1*eta1*(Theta'*Theta));
        F = mu1*Theta'*Y1;
        for i = 1:max_iters
            %X1 = X - mu1*(beta2*(B2*B2')*X + eta1*Theta'*(Theta*X - Y1));
            X1 = A*X+F;
            xn = norm(X - X1, 'fro');
            X = X1;
            if xn <= th
                if verbose
                    disp(['Converged in ' num2str(i) ' iterations']);
                end
                break;
            end
        end
    end
    
    out.X1 = X1;
end
