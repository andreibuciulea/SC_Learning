function w2_hat = check_smooth_X1(X1,B2,N)

    %%%%%%%%%%%%%%%%%%%%%%% Added by AGM
    %The purpose of this section is to identify the triangles which are active in our graph
    %Output: vector idxB2_true whose length is the number of triangles in
    %B2
    %First we generate the incididence matrices of a complete graph
    Boutxx = gen_B12(N); 
    B01xx = Boutxx.B1;B02xx = Boutxx.B2;
    %Secondly, we identify the indexes of the complety graph that correspond to our incidence matrix
    NTriang=size(B2,2);
    idxB2_true=zeros(1,NTriang);
    for t=1:NTriang
        idxxxx = (B02xx == B2(:,t));
        idxB2_true(t)=find(all(idxxxx));
    end
    
    
    %%AGM: Let's check how smooth signals are
    %We find the smoothness per triangle
    trian_smooth=sum((B02xx'*X1).^2,2);
    %To eliminate the triangles that do not involve edges that exist we set
    %the values of triangules 
    trian_smooth2=trian_smooth;
    trian_smooth2(trian_smooth<=1e-10)=1e6;
    %We now sort the triangules based on smoothness
    [trismo2_sort,idxtrismo2_sort]=sort(trian_smooth2','ascend');
    %Finally, we plot the results to check if the triangles with the smoothest signals are indeed the ones that are filled
    %Logical condition
    %disp(['True is the smoothest signals are the ones in filled triangles: ', num2str(all(sort(idxtrismo2_sort(1:NTriang))==sort(idxB2_true)))])
    figure
    stem(trian_smooth)
    hold on
    plot(idxB2_true,trian_smooth(idxB2_true),'rx')
    plot(idxtrismo2_sort(1:NTriang),trismo2_sort(1:NTriang),'gs')
    legend('Smoothness','Filled triangles','Triangles with smoothest signals')
    xlabel('Triangle index')
    ylabel('Level of smoothness')
    %axis([1 size(B2,2) 0.5 1.2*max(trian_smooth)])

    idxB2_hat = idxtrismo2_sort(1:NTriang);
    w2_true = zeros(size(B02xx,2),1);
    w2_hat = zeros(size(w2_true));
    w2_true(idxB2_true) = 1;
    w2_hat(idxB2_hat) = 1;
    w2_hat = logical(w2_hat);

    fi_tr = sum(w2_true);
    sm_tr_tr = sum((w2_true+w2_hat)==2);
    %disp(['-- From a total of ' num2str(fi_tr) ' true filled triangles ' num2str(sm_tr_tr) ' are the smoothest.' ])

%     figure()
%     plot(w2_true)
%     hold on
%     plot(w2_hat);




end