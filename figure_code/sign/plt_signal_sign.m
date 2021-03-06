% plots signal mean and SEM for the variability across sessions for
% sign/info/bursty

clear all 
close all
%clc 

savefig=0;
K=500;

%%

ntype={'sign','info','bursty'};
namebeh={'different','same'};

nameg={'minus','plus'};
namei=nameg;

task='plot signal sign sessions and regular/permuted';
disp(task);

savefile='/home/veronika/Dropbox/transfer/figure/sign/';
figname=['signal_sign_',sprintf('%1.0i',K)];

fs=10; % figure settings
lw=1.5;

blue=[0,0.48,0.74];
green=[0.2,0.6,0];
gray=[0.2,0.2,0.2];
red=[0.85,0.32,0.1];
col={blue,green};    

factor=100; % scaling for easier reading
pos_vec=[0,0,14,12];

%% load results

addpath(['/home/veronika/transfer_learning/result/signal/sign/K',sprintf('%1.0i',K)])

loadname=['sign_remove_info_',sprintf('%1.0i',K)];
load(loadname)

loadname2=['sign_perm_',sprintf('%1.0i',K)];
load(loadname2)

%%

x_mean=cellfun(@(x) squeeze(nanmean(x,2)).*factor, x_sign, 'UniformOutput', false);
x_diff=cellfun(@(x) squeeze(x(:,1,:)), x_mean, 'UniformOutput', false); % choice "different"
x_same=cellfun(@(x) squeeze(x(:,2,:)), x_mean, 'UniformOutput', false); % choice "same"

%% mean and standard error of the mean

K=size(x_same{1},2);                                                                     % number of time steps
nbses=size(x_same{1,1},1);                                                                      % number of sessions

mnm=cellfun(@(x) nanmean(x),x_diff,'UniformOutput',false);
mm=cellfun(@(x) nanmean(x),x_same,'UniformOutput',false);
snm=cellfun(@(x) nanstd(x),x_diff,'UniformOutput',false);
sm=cellfun(@(x) nanstd(x),x_same,'UniformOutput',false);

%% regular/permuted

regular=cellfun(@(x) squeeze(nanmean(x(:,2,:)- x(:,1,:))), x_mean, 'UniformOutput', false);
%permuted=cellfun(@(x) squeeze(nanmean(x(:,:,2,:)-x(:,:,1,:))).*factor, xperm, 'UniformOutput', false); % difference and average across sessions
if K==500
    permuted=cellfun(@(x) squeeze(nanmean(x(:,:,2,:)-x(:,:,1,:))).*factor, xp_sign, 'UniformOutput', false); % difference and average across sessions
else
    permuted=cellfun(@(x) squeeze(nanmean(x(:,:,2,:)-x(:,:,1,:))).*factor, xperm, 'UniformOutput', false); % difference and average across sessions
end
nperm=size(permuted{1},1);

%%

alpha=0.025/K;
idx=max([floor(alpha*nperm),1]);
lb=zeros(2,K);
ub=zeros(2,K);

for r=1:2
    for k=1:K
        x_sorted=sort(permuted{r}(:,k));
        lb(r,k)=x_sorted(idx);
        ub(r,k)=x_sorted(end-idx);
    end
    
end

%% plot

namesgn={'negative weight','positive weight'};

max_val=cat(1,cellfun(@(x,y) max(x+(y./nbses)),mm,sm),cellfun(@(x,y) max(x+(y./nbses)),mnm,snm));
maxy=max(max_val(:)).*2.5;

yt=-0.5:0.5:0.5;

xt=0:200:400;
x_vec=1:K;

H=figure('name',figname,'visible','on');

for r=1:2
    subplot(2,2,r)
    hold on
    
    % errorbars for the mean and the standard error of the mean(std(x)/sqrt(N))
    y1=mnm{r}-snm{r}./sqrt(nbses);
    y2=mnm{r}+snm{r}./sqrt(nbses);
    patch([x_vec fliplr(x_vec)], [y1 fliplr(y2)], col{1},'FaceAlpha',0.3,'EdgeColor',col{1})
    
    z1=mm{r}-sm{r}./sqrt(nbses);
    z2=mm{r}+sm{r}./sqrt(nbses);
    patch([x_vec fliplr(x_vec)], [z1 fliplr(z2)], col{2},'FaceAlpha',0.3,'EdgeColor',col{2})
    
    plot(x_vec,zeros(K,1),'--','color',gray)
    
    hold off
    box off
    
    ylim([-maxy,maxy])
    xlim([-2,K+2])
    
    set(gca,'YTick',yt)
    set(gca,'XTick',xt)
    
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    
    if r==1
        set(gca,'YTickLabel',yt)
    end
    
    title(namesgn{r},'FontName','Arial','fontsize',fs,'fontweight','normal')
    
    if r==1
        text(0.14,0.78,namebeh{1},'units','normalized','color',col{1},'FontName','Arial','fontsize',fs)
        text(0.14,0.9,namebeh{2},'units','normalized','color',col{2},'FontName','Arial','fontsize',fs)
        ylabel ('pop. signal (a.u.)','FontName','Arial','fontsize',fs);
    end
    
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
    
end

for r=1:2
    
    y1=lb(r,:);
    y2=ub(r,:);
    
    subplot(2,2,r+2)
    hold on
    patch([x_vec fliplr(x_vec)], [y1 fliplr(y2)], gray,'FaceAlpha',0.3,'EdgeColor',[0.7,0.7,0.7])
    plot(regular{r},'color',red,'linewidth',lw)
    plot(x_vec,zeros(K,1),'--','color',gray)
    
    hold off
    box off
    
    ylim([-maxy,maxy])
    xlim([-2,K+2])
    
    set(gca,'YTick',yt)
    set(gca,'XTick',xt)
    
    if r==1
        set(gca,'YTickLabel',yt, 'FontName','Arial','fontsize',fs)
    else
        set(gca,'YTickLabel',[])
    end
    
    if r==1
        text(0.14,0.9,'regular','units','normalized','color',red,'FontName','Arial','fontsize',fs)
        text(0.14,0.78,'permuted','units','normalized','color',gray,'FontName','Arial','fontsize',fs)
        ylabel ('difference of signals (a.u.)','FontName','Arial','fontsize',fs);
    end
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
end

axes
h1 = xlabel ('time (ms)','units','normalized','Position',[0.5,-0.08,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h1,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec)                                                      % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])      % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end
