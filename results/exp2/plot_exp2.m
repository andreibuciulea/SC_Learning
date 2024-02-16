load('results_exp2.mat')
merr = squeeze(median(err));
mrkrs = {'-s','-*','-^',':s',':*',':^'};
metrics = {'$L_{0}$','$L_{U}$'};
Models = {'SepSCL','RC','GreedySCL'};
nM = numel(Models);
plot_error = zeros(2,nM,numel(Ps));
%error values from Baseline and our model
plot_error(:,[1,3],:) = merr([3,4],:,:);

%error values for RC
%load('results_exp2_RC.mat','err');
load('res_RC.mat');
rc_error = squeeze(mean(err([3,5],:,:),2));
plot_error(:,2,:) = rc_error;
plot_error = plot_error/2;

figure('Position',[100,100,750,700])
i = 1;
for k = [1,2]
    for m = 1:nM
        lgd{i} = [metrics{k} '-' Models{m}];
        semilogy(Ps*100,squeeze(plot_error(k,m,:)),mrkrs{i},LineWidth=3,MarkerSize=14)
        hold on
        i = i+1;
    end
end
lg = legend(lgd,'FontSize',16,'FontWeight','bold', 'Interpreter','latex');
set(lg,'color','none');
xlabel('Percentage of available edge signals','FontSize',24,'FontWeight','bold', 'Interpreter','latex')
ylabel('Nerr for $L_{0}$ and $L_{U}$','FontSize',24,'FontWeight','bold', 'Interpreter','latex')
%title('Estimation error when varying the % of available edge signals')
xlim([55,95])
ylim([8e-3,1])
%ylim([1e-1,1e0])
yticks([1e-2,1e-1,1e0])
yticklabels({'10^{-2}','10^{-1}','10^{0}'})
xticks([55,60,65,70,75,80,85,90,95])
xticklabels({'55','60','65','70','75','80','85','90','95'})

ax = gca;
ax.FontSize = 19;
grid on
