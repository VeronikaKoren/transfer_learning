%% plot auc scores of single neurons in target and test 
% for V1 and V4

clear all
close all
format long

savefig=0;
info_case=1;

ba=1;
K=500;                                                           

%%
namea={'V1','V4'};
namei={'s+c','c'};

% figure settings
savefile='/home/veronika/Dropbox/transfer/figure/sign/';
figname=['weights_',namei{info_case},'_',sprintf('%1.0i',K)];

fs=11.5;
ms=5;
lw=1.2;
lwa=1;

gray=[0.7,0.7,0.7];
blue=[0,0.48,0.74];
red=[0.85,0.32,0.1];

baseline=0; 

%% load results

addpath('/home/veronika/transfer_learning/result/weight/weight_alltr/')

loadname=['weight_alltr_', namei{info_case},'_',namea{ba},'_',sprintf('%1.0i',K)];
load(loadname);

w=cell2mat(w_alltr);

%% nb of neurons with pos and neg weights

Ntot=length(w);

idx_pos=find(w > baseline);
npos=length(idx_pos);
mean_pos=mean(w(idx_pos));

idx_neg=find(w < baseline);
nneg=length(idx_neg);
mean_neg=mean(w(idx_neg));

display([npos,nneg]./Ntot, 'percent neur. with positive/negative weights')

%% test if weights are imbalanced (t-test)

[~,pval]=ttest(w);
display(pval,'pvalue ttest distribution of weights is different from 0')

%% plot 1 dist

col={blue,red};
[f,xi]=ksdensity(w,'Function','pdf');
fnorm=f./sum(f);

xt=[-1,0,1];
yt=[0,0.02];

pos_vec=[0,0,7,6];

H=figure('name',figname);
hold on
area(xi(1:50),fnorm(1:50),'FaceColor',col{1},'FaceAlpha',0.6)
area(xi(50:100),fnorm(50:100),'FaceColor',col{2},'FaceAlpha',0.6)
hold off

xlabel('weights','FontName','Arial','fontsize',fs)
ylabel('probability distribution','FontName','Arial','fontsize',fs)

axis([-1.4,1.4,0,0.033])
set(gca,'XTick',xt)
set(gca,'XTickLabel',xt,'FontName','Arial','fontsize',fs)
set(gca,'YTick',yt)
set(gca,'YTickLabel',yt,'FontName','Arial','fontsize',fs)

set(gca,'LineWidth',lwa,'TickLength',[0.025 0.025]);

set(H, 'Units','centimeters', 'Position', pos_vec) 
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)]) 

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

