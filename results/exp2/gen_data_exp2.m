clear all
addpath("..\opt")
addpath("..\utils")
load('RC_results_exp2_sparse.mat')
load('corr_data_exp2.mat');

[N,~,nG,ns]=size(all_H1);
outB12 = gen_B12(N);
B01 = outB12.B1;
B02 = outB12.B2;

figure()
subplot(221)
imagesc(all_H1(:,:,2,1))
subplot(222)
imagesc(all_C_X0_X1(:,:,2,1))
subplot(223)
imagesc(all_H1(:,:,1,9))
subplot(224)
imagesc(all_C_X0_X1(:,:,1,9))

figure()
subplot(221)
imagesc(all_H2(:,:,4,1))
subplot(222)
imagesc(all_H2(:,:,3,1))
subplot(223)
imagesc(all_H2(:,:,2,9))
subplot(224)
imagesc(all_H2(:,:,1,9))


%
err = zeros(5,nG,ns);
for g = 1:nG
    L0 = all_L0(:,:,g);
    L0d = L0-diag(diag(L0));
    A = diag(diag(L0))-L0;
    nA = norm(A,'fro')^2;
    nL0 = norm(L0d,'fro')^2;
   
    E = sum(all_w1(:,g));

    %compute laplacian of the true filled triangles
    T = sum(all_w2f(:,g));
    B2 = B02*diag(all_w2f(:,g));
    LB2 = B2*B2';
    LB2 = LB2-diag(diag(LB2));
    nLB2 =  norm(LB2,'fro')^2;


    for s = 1:ns

        %%%%% Error for the graph
        A_hat = all_H1(:,:,g,s);
        %measure error considering the obtained estimation
        err(1,g,s) = norm(A-A_hat,'fro')^2/nA;
        %measure error considering E highest links
        wm1 = all_wm1(:,g,s);
        Lm0 = B01*diag(wm1)*B01';
        Am0 = diag(diag(Lm0))-Lm0;
        A_hat(logical(Am0)) = 1;
        err(2,g,s) = norm(A-A_hat,'fro')^2/nA;
        L0_hat = diag(sum(A_hat))-A_hat;
        L0d_hat = L0_hat-diag(diag(L0_hat));
        err(3,g,s) = norm(L0d-L0d_hat,'fro')^2/nL0;
        [a,~] = sort(A_hat(:),'descend');
        A_hat_bin = A_hat >= a(E*2);
        err(4,g,s) = norm(A-A_hat_bin,'fro')^2/nA;
        
        %find w1_hat (edges of the estimated graph)
        w1_hat = compute_B1(A_hat_bin).w1;

%         figure()
%         subplot(221)
%         imagesc(A)
%         title(['A ' num2str(sum())])
%         colorbar()
%         subplot(222)
%         imagesc(A_hat)
%         title('A hat')
%         colorbar()
%         subplot(223)
%         imagesc(A_hat_bin)
%         title('A hat bin')
%         colorbar()
%         subplot(224)
%         imagesc(L0d_hat)
%         title('L0 hat')
%         colorbar()
        %%%% Error for the simplicial complexes
        %Find the B2 matrix of the triangles in T_hat 
        T_hat = all_T2_hat(:,:,g,s);
        T_vals = all_T2_vals(:,g,s);
        B2_hat = find_B2(T_hat,T_vals);

        %Using only the triangles associated with the estimated edges
        B2_aux = B2_hat;
        B2_aux(B2_aux>0) = 1;
        B2_aux(B2_aux<0) = -1;
        B2_aux = diag(w1_hat)*B2_aux;
        tr_idx = sum(abs(B2_aux))==3;
        B2_bin = B02*diag(tr_idx);
        B2_hat = B2_bin.*B2_hat;


        %select T triangles from B2_hat
        [~,b] = sort(sum(abs(B2_hat)),'descend');
        B2_hat_bin = B02(:,b(1:T));

        %only consider triangles between the estimated edges

        LB2_hat = B2_hat_bin*B2_hat_bin';
        LB2_hat = LB2_hat-diag(diag(LB2_hat));
        %compute error between B2_hat and 
        err(5,g,s) = norm(LB2-LB2_hat,'fro')^2/nLB2;

    end
end



