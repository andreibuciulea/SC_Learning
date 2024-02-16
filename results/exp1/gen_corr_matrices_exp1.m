addpath("..\..\opt")
addpath("..\.\utils")
nG = 100;
N = 20;
E0 = N*(N-1)/2;
T0 = trace((ones(N)-eye(N))^3)/6;
p = 0.4;%try to vary this
M0 = 1e3;%try to vary this
M1 = 50;
sig_type = 'smooth';
prms=struct('M0',M0,'M1',M1,'sig_type',sig_type);%parameters

all_L0 = zeros(N,N,nG);

%%% Full incidence matrices B01 and B02 
outB12 = gen_B12(N);
B01 = outB12.B1;B02 = outB12.B2; 

%%%for experiment 1
Sigmas = 0:0.03:0.5;
nS = numel(Sigmas);
all_C_X0 = zeros(N,N,nG,nS);
all_X0 = zeros(N,M0,nG,nS);


%%%for experiment 2
Ps = 0.8; %assuming partial observability of the edge signals and noisy node signals
nps = numel(Ps);
all_X1 = zeros(E0,M1,nG); % signals associated with all edges
all_Xm1 = zeros(E0,M1,nG,nps);% signals asscociated with selected edges
all_w1 = zeros(E0,nG);%all true edges
all_w2 = zeros(T0,nG);%all true triangles
all_w2f = zeros(T0,nG);%selected triangles
all_wm1 = zeros(E0,nG,nps);%edges associated with the observed edge signals

vars = zeros(5,nG,nS);

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
    

    %generate smooth edge signals
    X1 = generate_graph_signals1(B2*B2',prms).X_0;
    
    wm1 = set_M1(w1,Ps);%selection vector of the observed edge signals    

    %generate smooth node signals and assume corr = 1 for the observed edges
    Lm0 = B01*diag(wm1)*B01';
    Am0 = diag(diag(Lm0))-Lm0;

    for s = 1:nS
        %generate smooth and noisy node signals
        prms.sigma = Sigmas(s);
        out = generate_graph_signals0(L0,prms);
        X0 = out.X;
        all_X0(:,:,g,s) = X0;
        C_X0 = corr(X0');

        vars(1,g,s) = out.snr;
        vars(2,g,s) = out.px;
        vars(3,g,s) = out.pn;
        vars(4,g,s) = out.sigma_n;
        

        %C_X0(logical(Am0)) = 0; %for Rc these entries should be zero to consider them as observed edges
        all_C_X0(:,:,g,s) = C_X0;
    end

    %store all the necesary data
    all_L0(:,:,g) = L0; 
    all_w1(:,g) = w1;
    all_wm1(:,g) = wm1;
    all_w2(:,g) = w2;
    all_w2f(:,g) = w2_f;


    all_X1(:,:,g) = diag(w1)*(X1);
    all_Xm1(:,:,g) = diag(wm1)*X1; %observed edge signals

end

%save('corr_data_exp1.mat','all_C_X0','all_L0','all_w1','all_w2','all_w2f','all_wm1','all_X0','all_X1','all_Xm1','Sigmas');