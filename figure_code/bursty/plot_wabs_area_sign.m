%% test if neurons with strong weight have small area under the PS (a.k.a. bursty neurons have stronger weights)

close all
clear all

savefig=0;

%%

ba=1;
period=2;
K=400;

namei={'s+c','c'};
namea={'V1','V4'};
namep={'tar','test'};
figname=['wabs_area_sign_',sprintf('%1.0i',K)];
savefile='/home/veronika/Dropbox/transfer/figure/bursty/';
pos_vec=[0,0,7,12];

fs=11;
ms=7;
lw=1.2;
lwa=1;

red=[0.85,0.32,0.1];
blue=[0,0.48,0.74];

col={red,blue};

%% load area under PS

addpath('/home/veronika/synced/transfer_result/basic_stat/ps/')
loadname=['ps_',namea{ba},namep{period},'_',sprintf('%1.0i',K)];
load(loadname)

loadname2=['ps_poisson_',namea{ba},namep{period},'_',sprintf('%1.0i',K)];
load(loadname2)

vec=3:41;
areaps=cell2mat(cellfun(@(x) nanmean(x(:,vec),2), ps_sess,'UniformOutput', false)); % area regular
areaps0=cell2mat(cellfun(@(x) squeeze(nanmean(x(:,:,vec),3))', psp,'UniformOutput', false)); % area permuted

%% load weights

addpath('/home/veronika/synced/transfer_result/weight/weight_alltr/')

loadname=['weight_alltr_s+c','_V1_', sprintf('%1.0i',K)];
load(loadname)

loadname0=['weight_allp_s+c','_V1_', sprintf('%1.0i',K)];
load(loadname0)

wabs=abs(cell2mat(w_alltr)); % strength of the weight regular
wabs0=cell2mat(cellfun(@(x) abs(x)',w_allp,'UniformOutput', false)); % strength of the weight permuted

%% compute the correlation coefficient and the p-value

tags=[1,-1];

pval=zeros(2,1);
Rround=zeros(2,1);
for s=1:2
    
    wsgn=sign(cell2mat(w_alltr))== tags(s);
    idx_use=find(wsgn);
    
    % round to 3 decimals
    
    R=corr(wabs(idx_use),areaps(idx_use),'tail','left','type','Pearson');
    Rround(s)=round(R*1000)/1000;
   
    % R of permuted
    
    nperm=size(wabs0,2);
    %N=size(wabs0,1);
    R0=zeros(nperm,1);
    
    for p=1:nperm
        R0(p)=corr(wabs0(idx_use,p),areaps0(idx_use,p));
    end
    
    pval(s)=sum(R0<R)/nperm;
    
end

display(Rround,'linear correlation R(abs(w),area_ps)')
display(pval,'p-value permutation test')

%% plot

titles={'positive weight','negative weight'};

xt=0:0.5:1;
yt=0.8:0.2:1.2;

H=figure('name',figname,'visible','on');

for s=1:2
    
    wsgn=sign(cell2mat(w_alltr))== tags(s);
    idx_use=find(wsgn);
    N_group=length(idx_use)
    
    subplot(2,1,s)
    
    scatter(wabs(idx_use),areaps(idx_use),10,'MarkerFaceColor',col{s},'MarkerEdgeColor',col{s})
    line([0,1],[1,1],'color',[0.5,0.5,0.5])
    ls=lsline;
    ls.Color='k';
    if sign(Rround(s))==1
        text(0.5,0.8,['R = ',sprintf('%0.3f',Rround(s))],'units','normalized', 'FontName','Arial','Fontsize',fs)
    else
        text(0.5,0.8,['R = - ',sprintf('%0.3f',abs(Rround(s)))],'units','normalized', 'FontName','Arial','Fontsize',fs)
    end
    if pval(s) < 0.05/2
        text(0.45,0.8,'*','units','normalized','color','red','FontName','Arial','Fontsize',fs+5)
    end
    grid on
    
    axis([0,1,0.7,1.3])
    title(titles{s},'FontName','Arial','Fontsize',fs,'Fontweight','normal')
    
    set(gca,'YTick',yt)
    set(gca,'YTickLabel',yt, 'FontName','Arial','Fontsize',fs)
    set(gca,'XTick',xt)
    
    if s==2
        xlabel('strength weight','FontName','Arial','Fontsize',fs)
        set(gca,'XTickLabel',xt,'FontName','Arial','Fontsize',fs)
    else
        set(gca,'XTickLabel',[])
    end
    
    op=get(gca,'OuterPosition');
    set(gca,'OuterPosition',[op(1)+0.06 op(2) op(3)-0.06 op(4)]); % OuterPosition = [left bottom width height]
    set(gca,'LineWidth',lwa,'TickLength',[0.025 0.025]);
    
end

axes
h2 = ylabel ('area under power spectrum','units','normalized','Position',[-0.07,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h2,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec)
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)]) % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end
