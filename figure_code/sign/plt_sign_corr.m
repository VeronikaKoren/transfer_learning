% plots cross-correlation function for subnetworks

clear all 
close all
clc 

savefig=0;
K=400;

namei={'minus','plus'};
variable={'r_sign','r_info','r_ps'};

task='plot sign correlation';
disp(task);

savefile='/home/veronika/Dropbox/transfer/figure/sign/';
figname=['sign_corr_',sprintf('%1.0i',K)];

pos_vec=[0,0,9,12];
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

%% load results

addpath(['/home/veronika/transfer_learning/result/signal/sign/K',sprintf('%1.0i',K)])
xcorr_tau=cell(3,1);
for i=1:3
    loadname=['xcorr_sign_tau', sprintf('%1.0i',tau_range(i))]; % correlation function in sessions, both conditions together
    load(loadname);
    xcorr_tau{i}=r_sign;
    
end

%% mean and std across sessions

L=size(xcorr_tau{1},2);
K=(L+1)/2;
lags=-K+1:K-1;
nbses=size(xcorr_tau{1},1);

mm=cellfun(@(x) nanmean(x),xcorr_tau,'UniformOutput',false);
sm=cellfun(@(x) nanstd(x),xcorr_tau,'UniformOutput',false);

%% permuted tau=20;

loadname2=['xcorr_perm_tau', sprintf('%1.0i',tau_range(2))]; % correlation function in sessions, both conditions together
load(loadname2);

nperm=size(r_perm,2);
rp=squeeze(nanmean(r_perm)); % mean across sessions
%%

alpha=0.025/K;
idx=max([floor(alpha*nperm),1]);
        
lb=zeros(1,L);
ub=zeros(1,L);
for l=1:L
    x=sort(squeeze(rp(:,l)));
    lb(l)=x(idx);
    ub(l)=x(end-(idx-1));
end

%% plot

mini=-0.05;
maxi=0.01;
yt=-0.04:0.02:0;
xt=-300:300:300;
ytext=0.34;

pltidx=[3,2,1];
nlag=length(lags);

H=figure('name',figname);
subplot(2,1,1)
for i=1:3
    
    hold on
    
    z1=mm{pltidx(i)}-sm{pltidx(i)}./sqrt(nbses);
    z2=mm{pltidx(i)}+sm{pltidx(i)}./sqrt(nbses);
    patch([lags fliplr(lags)], [z1 fliplr(z2)], col{pltidx(i)},'FaceAlpha',0.3,'EdgeColor',col{pltidx(i)})
    
    plot(lags,zeros(nlag,1),'--','color',gray)
    
    hold off
    box off
    
    text(0.15,ytext + 0.12*(i-1),['\lambda=',namelambda{i}],'color',col{i},'units','normalized', 'FontName','Arial','fontsize',fs)
    
    ylim([mini,maxi])   
    xlim([-K,K])
    
    
end
set(gca,'YTickLabel',yt)
grid on
set(gca,'YTick',yt)
set(gca,'YTickLabel',yt, 'fontname','Arial', 'fontsize',fs)
set(gca,'XTick',xt)    
set(gca,'XTickLabel',[])
set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
op=get(gca,'OuterPosition');
set(gca,'OuterPosition',[op(1)+0.06 op(2) op(3)-0.05 op(4)]); % OuterPosition = [left bottom width height]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,2)

hold on
patch([lags fliplr(lags)], [lb fliplr(ub)], gray,'FaceAlpha',0.5,'EdgeColor',gray)
plot(lags,zeros(nlag,1),'--','color',gray)
plot(lags,mm{2},'color',col{2},'linewidth',lw)

ylim([mini,maxi])
xlim([-K,K])

hold off
box off

text(0.15,ytext,'regular','color',col{2},'units','normalized', 'fontname','Arial', 'fontsize',fs)
text(0.15,ytext-0.13,'permuted','color',gray,'units','normalized', 'fontname','Arial', 'fontsize',fs)

op=get(gca,'OuterPosition');
set(gca,'OuterPosition',[op(1)+0.06 op(2) op(3)-0.05 op(4)]); % OuterPosition = [left bottom width height]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gca,'YTickLabel',yt, 'fontname','Arial', 'fontsize',fs)
grid on
set(gca,'YTick',yt)
set(gca,'XTick',xt)    
set(gca,'XTickLabel',xt, 'fontname','Arial', 'fontsize',fs)
set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes

h1 = xlabel ('time lag (ms)','units','normalized','Position',[0.5,-0.07,0],'FontName','Arial','fontsize',fs);
h2 = ylabel ('correlation function','units','normalized','Position',[-0.08,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h2,'visible','on')
set(h1,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec) % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%%
peak=cellfun(@(x) x(K),mm)



