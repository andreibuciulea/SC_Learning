function out = compute_B2(A,B,w1)
    [N,E] = size(B);
    T = trace(A^3)/6;
    B2 = zeros(E,T);
    w2 = zeros(T,1);
    option = 2;
    if option == 1
        NN = 1:N;
        t = 1;
        for e = 1:E
            [idx] = find(B(:,e)~=0);
            if ~isempty(idx)
                NN = 1:N;
                NN([idx(1),idx(2)]) = [];
                for n = NN
                    if A(idx(1),n)*A(idx(2),n) && (n>idx(2))
                        idx(3) = n;
                        idx = sort(idx);
                        idxL(1) = find(B(idx(1),:).*B(idx(2),:)~=0); 
                        idxL(2) = find(B(idx(1),:).*B(idx(3),:)~=0); 
                        idxL(3) = find(B(idx(2),:).*B(idx(3),:)~=0); 
                        B2(idxL([1,3]),t) = 1;
                        B2(idxL(2),t) = -1;
                        t = t+1;
                        %disp(['Nodes ' num2str(idx(1)) num2str(idx(2)) num2str(idx(3)) ' form a triangle with link idx: [' num2str(idxL(1)) ',' num2str(idxL(2)) ',' num2str(idxL(3)) ']'])
                    end
                end
            end
        end
        B2=flip(unique(B2','rows')',2);
    elseif option == 2
        [~,w1_idx] = find(w1~=0); 
        B02 = gen_B12(N).B2;
        B2_sel = diag(w1)*B02;
        w2 = sum(abs(B2_sel))==3; %binary vector with 1 for triangles
        [~,w2_idx] = find(sum(abs(B2_sel))==3);%index of the triangles
        B2 = B02(:,w2_idx);
        out.B2c = diag(w1)*B02*diag(w2);
    end
    %disp(t)
    out.B2 = B2;
    out.w2 = w2';
end