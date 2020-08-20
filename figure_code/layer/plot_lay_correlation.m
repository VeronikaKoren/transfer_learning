% plots cross-correlation function of + and - neurons in V1 and V4 during tar, delay and test

clear all 
close all
clc 

savefig=0;

task='plot layer correlation';
savefile='/home/veronika/Dropbox/transfer/figure/layer/';
figname='lay_corr';

pos_vec=[0,0,14,16];
fs=12; 
lw=1.5;

blue=[0,0.48,0.74];
red=[0.85,0.32,0.1];
gray=[0.2,0.2,0.2];

col={'m',red,blue};

tau_range=[10,20,40];
namelambda=cell(length(tau_range),1);
for t=1:3
    namelambda{t}=[sprintf('%1.0i',tau_range(t)),'^{-1}'];
end

%% load data

addpath('/home/veronika/transfer_learning/result/signal/layer/')

% regular model for different time constant lambda
xcorr_tau=cell(3,3);

for t=1:3

    loadname=['xcorr_layer_tau', sprintf('%1.0i',tau_range(t))]; % correlation function in sessions, both conditions together
    load(loadname);
    for i=1:3
        xcorr_tau{i,t}=squeeze(r_layer(:,i,:));
    end
    
end

%% mean and std across sessions

L=size(xcorr_tau{1,1},2);
K=(L+1)/2;
lags=-K+1:K-1;
nbses=size(xcorr_tau{1,1},1);

mm=cellfun(@(x) nanmean(x),xcorr_tau,'UniformOutput',false);
sm=cellfun(@(x) nanstd(x),xcorr_tau,'UniformOutput',false);

%% permuted tau=20;

loadname2=['xcorr_layer_perm_tau', sprintf('%1.0i',tau_range(2))]; % correlation function in sessions, both conditions together
load(loadname2);
nperm=size(r_layer_perm,3);

rp=squeeze(nanmean(r_layer_perm)); % upper and lower bound

%%

alpha=0.025/K;
idx=max([floor(alpha*nperm),1]);
        
lb=zeros(3,L);
ub=zeros(3,L);
for t=1:3
    for l=1:L
        x=sort(squeeze(rp(t,:,l)));
        lb(t,l)=x(idx);
        ub(t,l)=x(end-(idx-1));
    end
end

%% plot

mini=-0.01;
maxi=0.035;
yt=0:0.02:0.02;
xt=-300:300:300;

pltidx=[1,3,5];

H=figure('name',figname);
for i=1:3
    
    subplot(3,2,pltidx(i))
    hold on
    
    for t=1:3
        
        z1=mm{i,t}-sm{i,t}./sqrt(nbses);
        z2=mm{i,t}+sm{i,t}./sqrt(nbses);
        patch([lags fliplr(lags)], [z1 fliplr(z2)], col{t},'FaceAlpha',0.3,'EdgeColor',col{t})
        
        plot(lags,zeros(L,1),'--','color',gray)
        box off
        
        ylim([mini,maxi])
        xlim([-K,K])
        
        if i==1
            text(0.05,0.6 + 0.12*(t-1),['\lambda=',namelambda{t}],'color',col{t},'units','normalized', 'FontName','Arial','fontsize',fs)
        end
        
    end
    hold off
    grid on
    
    
    yticks(yt)
    xticks(xt)
    
    ax=gca;
    if i==3
        ax.XTickLabel=xt;
    else  
        ax.XTickLabel={};
    end
    set(gca,'YTickLabel',yt, 'FontName','Arial','fontsize',fs)
    
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
    subplot(3,2,pltidx(i)+1)
    
    plot(lags,mm{i,2},'color',col{2},'linewidth',lw)
    hold on
    patch([lags fliplr(lags)], [lb(i,:) fliplr(ub(i,:))], gray,'FaceAlpha',0.5,'EdgeColor',gray)
    
    ylim([mini,maxi])
    xlim([-K,K])
    
    hold off
    box off
    
    text(0.97,0.5,name2lay{i},'units','normalized', 'FontName','Arial','fontsize',fs,'rotation',25)
   
    if i==1
        text(0.1,0.85,'regular','color',col{2},'units','normalized','FontName','Arial','fontsize',fs)
        text(0.1,0.75,'permuted','color',gray,'units','normalized', 'FontName','Arial','fontsize',fs)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold off
    grid on
    
    yticks(yt)
    xticks(xt)
    
    ax=gca;
    set(gca,'YTickLabel', yt, 'FontName','Arial','fontsize',fs)
    ax.YTickLabel={};
    if i==3
        ax.XTickLabel=xt;
    else  
        ax.XTickLabel={};
    end
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
    
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes
%text(-0.15,1.05,'B','units','normalized','FontName','Arial','fontsize',fs,'fontweight','bold')
h1 = xlabel ('time lag (ms)','units','normalized','Position',[0.5,-0.07,0],'FontName','Arial','fontsize',fs);
h2 = ylabel ('correlation function','units','normalized','Position',[-0.1,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h2,'visible','on')
set(h1,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec) % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%%
peak=cellfun(@(x) x(K),mm);
display(peak)


