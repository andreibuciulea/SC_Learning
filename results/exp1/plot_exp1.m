load('results_ours_v2.mat')
merr = squeeze(mean(err));
mrkrs = {'-s','-*','-^',':s',':*',':^'};
metrics = {'$L_{0}$','$L_{U}$'};
Models = {'SepSCL','RC','GreedySCL'};
nM = numel(Models);
plot_error = zeros(2,nM,numel(Sigmas));
%error values from Baseline and our model
plot_error(:,[1,3],:) = merr([3,4],:,:);

%error values for RC
%load('results_RC.mat','err');
load('res_RC.mat');

rc_error = squeeze(mean(err([3,5],:,:),2));
plot_error(:,2,:) = rc_error;

plot_error = plot_error/2;
sv = 1:numel(Sigmas);% [1, 3, 5:3:26];
figure('Position',[100,100,750,700])
i = 1;
for k = [1,2]
    for m = 1:nM
        lgd{i} = [metrics{k} '-' Models{m}];
        semilogy(Sigmas(sv),squeeze(plot_error(k,m,sv)),mrkrs{i},LineWidth=3,MarkerSize=14)
        hold on
        i = i+1;
    end
end
lg = legend(lgd,'FontSize',16,'FontWeight','bold', 'Interpreter','latex');
set(lg,'color','none');
xlabel('Normalized noise error','FontSize',24,'FontWeight','bold', 'Interpreter','latex')
ylabel('Nerr for $L_{0}$ and $L_{U}$','FontSize',24,'FontWeight','bold', 'Interpreter','latex')
%title('Recovery error with 80% of the edges known')
ylim([1e-2,1])
xlim([0,0.3])
%xlim([0,max(msnr)]3
grid on

%ylim([1e-1,1e0])
%yticks([0.1,0.2,0.5,1,2])
%yticklabels({'0.1','0.2','0.5','1','2'})
%xticks([0,2,4,6,8,10,11])
%xticklabels({'0','2','4','6','8','10','11'})

ax = gca;
ax.FontSize = 19;
grid on


msnr = 10*log10(mean(snr(:,sv)));
msnr(1) = 12;
figure('Position',[100,100,850,600])
i = 1;
for k = [1,2]
    for m = 1:nM
        lgd{i} = [metrics{k} ' ' Models{m}];
        semilogy(msnr,squeeze(plot_error(k,m,sv)),mrkrs{i},LineWidth=3,MarkerSize=14)
        hold on
        i = i+1;
    end
end
lg = legend(lgd,'FontSize',16,'FontWeight','bold', 'Interpreter','latex');
set(lg,'color','none');
xlabel('Signal to noise ratio (dB)','FontSize',24,'FontWeight','bold', 'Interpreter','latex')
ylabel('Nerr for $L_{0}$ and $L_{U}$','FontSize',24,'FontWeight','bold', 'Interpreter','latex')
%title('Recovery error with 80% of the edges known')
ylim([0.03,2])
xlim([0,12])
%xlim([0,max(msnr)])
grid on

%ylim([1e-1,1e0])
%yticks([0.1,0.2,0.5,1,2])
%yticklabels({'0.1','0.2','0.5','1','2'})
xticks([0,2,4,6,8,10,12])
xticklabels({'0','2','4','6','8','10','12'})

ax = gca;
ax.FontSize = 19;
grid on

