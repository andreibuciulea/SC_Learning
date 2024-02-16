addpath("..\opt")
addpath("..\utils")
clear all
nG = 100;
N = 20;
E0 = N*(N-1)/2;
T0 = trace((ones(N)-eye(N))^3)/6;
p = 0.3;%try to vary this
M0 = 1e3;%try to vary this
M1 = 50;
sig_type = 'smooth';
prms=struct('M0',M0,'M1',M1,'sig_type',sig_type);%parameters

all_L0 = zeros(N,N,nG);

%%% Full incidence matrices B01 and B02 
outB12 = gen_B12(N);
B01 = outB12.B1;B02 = outB12.B2; 


%%%for experiment 2
Ps = 0.95:-0.05:0.55;
nps = numel(Ps);
all_C_X0_X1 = zeros(N,N,nG,nps);
all_X0 = zeros(N,M0,nG);
all_X1 = zeros(E0,M1,nG); % signals associated with all edges
all_Xm1 = zeros(E0,M1,nG,nps);% signals asscociated with selected edges
all_w1 = zeros(E0,nG);%all true edges
all_w2 = zeros(T0,nG);%all true triangles
all_w2f = zeros(T0,nG);%selected triangles
all_wm1 = zeros(E0,nG,nps);%edges associated with the observed edge signals


for g = 1:nG

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
    B2 = B02(:,logical(w2_f));%B2 of the selected triangles as filled
    

    %generate smooth node signals and correlation matrix
    X0 = generate_graph_signals0(L0,prms).X;
    C_X0 = corr(X0');

    %generate smooth edge signals
    X1 = generate_graph_signals1(B2*B2',prms).X_0;

    
    for r = 1:nps
        ps = Ps(r); %ratio of observed edge signals [0,1]
        wm1 = set_M1(w1,ps);%selection vector of the observed edge signals    
        all_Xm1(:,:,g,r) = diag(wm1)*X1; %observed edge signals
        all_wm1(:,g,r) = wm1; % all the selection vectors of the observed edge signals

        %generate smooth node signals and assume corr = 1 for the observed edges
        Lm0 = B01*diag(wm1)*B01';
        Am0 = diag(diag(Lm0))-Lm0;
        %C_X0(logical(Am0)) = 0.6;
        all_C_X0_X1(:,:,g,r) = C_X0;%all the covariances of X0 varying the number of observed edges
    end

    
    %store all the necesary data
    all_L0(:,:,g) = L0; %all Laplacian matrices
    all_w1(:,g) = w1; %all edge selection vectors
    all_w2(:,g) = w2; %all triangle selection vectors
    all_w2f(:,g) = w2_f;%all filled triangle selection vectors
    all_X0(:,:,g) = X0; %we generate one signal X for each graph and modify the number of observed edges
    all_X1(:,:,g) = diag(w1)*(X1);%all true edge signals
    

end
save('corr_data_exp2.mat','all_C_X0_X1','all_L0','all_w1','all_w2','all_w2f','all_wm1','all_X0','all_X1','all_Xm1','Ps');