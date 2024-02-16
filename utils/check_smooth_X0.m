function w1s = check_smooth_X0(X0,w1,N) 
    B01xx= gen_B12(N).B1;
    E0 = size(B01xx,2);
%%%%%%%%%%%%%%Check edge smoothness 
    C = X0*X0';
    c2 = diag(B01xx'*C*B01xx*eye(E0));
    smooth = c2;
    [a,~] = sort(c2);c2 = c2<=a(sum(w1));c2 = c2>=1e-5;
    w1s = c2;
    figure(1)
    subplot(211)
    stem(0.2*w1)
    hold on
    plot(smooth/max(max(smooth)))
    legend('True edges','Smoothness in each edge')
    subplot(212)
    stem(w1)
    hold on
    stem(c2)
    legend('True edges','Smoothest edges')

    t_s = sum((c2+w1)==2);
    t_e = sum(w1);
    %disp(['-- From a total of ' num2str(t_e) ' true edges ' num2str(t_s) ' are the smoothest.' ])
end