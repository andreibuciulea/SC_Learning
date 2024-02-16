addpath("opt")
addpath("utils")
clear all
close all

%clc
%rng(2)

N = 20; %number of nodes
E0 = N*(N-1)/2; %numer of all possible edges
M0 = 1e2; %number of node signals
M1 = 50; %number of edge signals
sig_type = 'smooth';%type of the signal
Models = {'check','ours',};%algoritm type for graph estimation
nM = numel(Models);

p = 0.4;%link probability for ER graphs

%%% Full incidence matrices B01 and B02 
outB12 = gen_B12(N);
B01 = outB12.B1;B02 = outB12.B2; 

tic

% params for estimate the graph topology
a0 = 1;a1 = 1;a2 = 1;
rho = 0.1;ga = 0.1;
la0 = 1e3;la1 = 1e3;
max_iters = 5;
verbose = false;
input = struct('a0',a0,'a1',a1,'a2',a2,'max_iters',max_iters,...
               'la0',la0,'la1',la1,'ga',ga,'rho',rho,'verbose',verbose);

nG = 32;
Sigmas = 0:0.05:0.5;
ns = numel(Sigmas);
ps = 0.80; %ratio of observed edge signals [0,1]
err = zeros(nG,5,nM,ns);
snr = zeros(nG,ns);
tic

parfor g = 1:nG
    g
    err_g = zeros(5,nM,ns);
        
    %generate graph
    A = generate_connected_ER(N,p);%adjacency matrix
    L0 = diag(sum(A))-A; %Laplacian matrix
    outB1 = compute_B1(A);
    w1 = outB1.w1; B1 = outB1.B1;%incidence matrix B1
    outB2 = compute_B2(A,B1,w1);
    w2 = outB2.w2; B2 = outB2.B2; %incidence matrix B2

    %select half of the triangles as filled
    T = trace(A^3)/6; %number of triangles
    w2_idx = find(w2==1);
    Tidx = sort(randperm(T,round(T/2)));
    w2_f = zeros(size(w2));
    w2_f(w2_idx(Tidx)) = 1;
    B2 = B02(:,logical(w2_f));

    if verbose
        disp("%%%%%%%%%%%%%%%%%%%%%")
        disp(['Number of edges: ' num2str(sum(w1))])
        disp(['Number of filled triangles: ' num2str(size(B2,2))])
    end
    L1 = B1'*B1+B2*B2'; %hodge Llaplacian
    LB1 = B1'*B1;
    LB2 = B2*B2';


    %generate smooth edge signals
    %params for generate smooth edge signals
    prms=struct('M0',M0,'M1',M1,'sig_type',sig_type);%parameters
    tmp = generate_graph_signals1(B2*B2',prms);
    X_B2 = tmp.X_0;
    Y_B2 = tmp.X;
    X1 = diag(w1)*(X_B2);
    Y1 = diag(w1)*(Y_B2);
    
    %Set a mask M1 for edge signals
    wm1 = set_M1(w1,ps);%selection vector of the observed edge signals    
    Xm1 = diag(wm1)*X1; %observed edge signals

    %find the mask for the unobserved edges
    Lm0 = B01*diag(wm1)*B01';
    Am0 = diag(diag(Lm0))-Lm0;
    nAm0 = not(logical(Am0));
    Ao = nAm0.*A;
    nAo = norm(Ao,'fro')^2;    
    nA = norm(A,'fro')^2;



    for r = 1:ns
        %generate smooth node signals
        %params for generate smooth graph signals
        prms=struct('M0',M0,'M1',M1,'sig_type',sig_type,'sigma',Sigmas(r));%parameters
        tmp = generate_graph_signals0(L0,prms);
        X0 = tmp.X_0;
        Y0 = tmp.X;
        snr(g,r) = tmp.snr;
        %Set a mask M0 for node signals
        %M0 = ones(size(X0));
        %Xm0 = M0.*Y0; 
        
        %%%add variables in input structure
        %input.w1 = w1;
        %input.B2 = B2;
    
        %estimate the graph and SCs from signals
        T = size(B2,2); %number of the filled triangles
        E = sum(w1); %number of the actual edges
        %input.B1 = B1;
        
        if verbose
            disp(['Percentage of observed edge signals: ' num2str(ps*100)])
            %check_smooth_X1(X1,B2,N);
            %check_smooth_X1(Xm1,B2,N);
            %check_smooth_X0(X0,w1,N);
        end

        for m = 1:nM 
            model = Models{m};

            out = estimate_SC(Y0,Xm1,E,T,model,input,B2,w1);
            
            %compute the estimation error
            B1_hat = out.B1; B2_hat = out.B2;X1_hat = out.X1;
            L0_hat = B1_hat*B1_hat'; LB1_hat = B1_hat'*B1_hat;LB2_hat = B2_hat*B2_hat';
            w1_hat = out.w1; w2_hat = out.w2;
            L1_hat = LB1_hat+LB2_hat;
            A_hat = diag(diag(L0_hat))-L0_hat;
            Ao_hat = nAm0.*(A_hat);
            err_g(1,m,r) = norm(B1-B1_hat,'fro')^2/norm(B1,'fro')^2;
            err_g(2,m,r) = norm(A-A_hat,'fro')^2/nA;
            err_g(3,m,r) = norm(L0-L0_hat,'fro')^2/norm(L0,'fro')^2;
            err_g(4,m,r) = norm(LB2-LB2_hat,'fro')^2/norm(LB2,'fro')^2;
            err_g(5,m,r) = norm(Ao-Ao_hat,'fro')^2/(nAo+norm(Ao_hat,'fro')^2);
            %err_g(5,m,r) = norm(X1-X1_hat,'fro')^2/norm(X1,'fro')^2;
        
            if verbose
                t_s = sum((w1+w1_hat)==2);
                disp(['++' model '++ From a total of ' num2str(E) ' true edges ' num2str(t_s) ' are the smoothest.' ])
            
                sm_tr_tr = sum((w2_f+w2_hat)==2);
                disp(['++' model '++ From a total of ' num2str(T) ' true filled triangles ' num2str(sm_tr_tr) ' are the smoothest.' ])
            end
        end
    end
    err(g,:,:,:) = err_g;
end
toc
% figure
% subplot(231)
% imagesc(L0)
% title('L0')
% colorbar()
% subplot(232)
% imagesc(LB1)
% title('LB1')
% colorbar()
% subplot(233)
% imagesc(LB2)
% title('LB2')
% colorbar()
% subplot(234)
% imagesc(L0_hat)
% title('L0 hat')
% colorbar()
% subplot(235)
% imagesc(LB1_hat)
% title('LB1 hat')
% colorbar()
% subplot(236)
% imagesc(LB2_hat)
% title('LB2 hat')
% colorbar()

%% plot results
merr = squeeze(mean(err));

mrkrs = {'-s',':^','--x','-s',':^','--x'};
metrics = {'B1','A-hat','L0','LB2','Ao hat'};
figure()
i = 1;
for k = [3,4]
    for m = 1:nM
        lgd{i} = [metrics{k} ' ' Models{m}];
        semilogy(Sigmas,squeeze(merr(k,m,:)),mrkrs{m},LineWidth=2)
        hold on
        i = i+1;
    end
end
legend(lgd)
xlabel('Normalized noise error')
ylabel('Normalized error')
grid on
