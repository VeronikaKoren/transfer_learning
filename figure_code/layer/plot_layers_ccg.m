
clear all
close all
clc

savefig=0;

ba=1;
period=2;
K=400;
%%

namea={'V1','V4'};
namep={'target','test'};
nameg={'SG','G','IG'};
ng=length(nameg);

figname=['ccg_layer_',sprintf('%1.0i',K)];
savefile='/home/veronika/Dropbox/transfer/figure/layer/ccg/';

red=[0.85,0.32,0.1];
blue=[0,0.48,0.74];
magenta=[1,0,1];
black=[0,0,0];
col={magenta,blue, black};

fs=9;
ms=5;
lw=1.2;
lwa=1;

%% plot ccg sign for plus and minus

addpath('/home/veronika/transfer_learning/result/pairwise/')
loadname=['ccg_layer_',sprintf('%1.0i',K),'.mat'];                
load(loadname)

%%

iW=(size(lags,2)+1)/2;
speak=cellfun(@(x) x(:,iW),ccg3,'un',0);
display(pval,'p-value groups');

xvec=-iW+1:iW-1;
msyn=cellfun(@mean,speak);                                    % average synchrony
semsyn=cellfun(@(x) std(x)./sqrt(length(x)),speak);

mc=cellfun(@(x) mean(x,2), r_all,'UniformOutput', false);    % average rccg acros neurons
semc=cellfun(@(x) std(x,0,2)./sqrt(size(x,2)),r_all,'UniformOutput', false); % std across neurons

for c=1:3
    strcat('pval_', nameg{idxes{1}(c)},' &_ ', nameg{idxes{2}(c)},' is_', sprintf('%0.4f', pval(c)))
 
end

%% plot

pos_vec=[0 0 18 5];
titles={'Cross-corr.','Synchrony','Average synchrony','Area under ccg'};

yt=0:0.04:0.08;
ylimit=[-0.02,0.08];
dx=abs(ylimit(2)-ylimit(1))/20;

H=figure('name',figname,'visible','on');

%%%%%%%%%%%%%%%%%%%   ccg with median peak for eack layer
xt=[-100,0,100];
subplot(1,4,1)   

hold on
for ii=[3,2,1]
    idx=find(speak{ii}==median(speak{ii}));
    plot(xvec,ccg3{ii}(idx,:),'color',[col{ii},0.5],'linewidth',lw+0.5);
    text(0.7,0.9-(ii-1)*0.1,nameg{ii},'units','normalized','FontName','Arial','fontsize',fs,'color',col{ii})
end
hold off

xlim([xt(1),xt(end)])
ylim(ylimit)
grid on

title(titles{1}, 'FontName','Arial','Fontsize',fs,'Fontweight','normal')

xlabel('lag (ms)', 'FontName','Arial','Fontsize',fs)
ylabel ('Prob. coincidence','FontName','Arial','fontsize',fs)

set(gca,'YTick',yt)
set(gca,'YTickLabel',yt, 'FontName','Arial','Fontsize',fs)
set(gca,'XTick',xt)
set(gca,'XTickLabel',xt,'FontName','Arial','Fontsize',fs)

set(gca,'LineWidth',lwa,'TickLength',[0.025 0.025]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eCDF of
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% synchrony

yt=0:0.5:1;

subplot(1,4,2)
hold on
for ii=1:ng
    [f,xi]=ecdf(speak{ii});
    plot(xi,f,'color',col{ii},'linewidth',lw)
end
hold off

grid on
xlim([-0.05,0.25])

title(titles{2}, 'FontName','Arial','Fontsize',fs,'Fontweight','normal')
xlabel ('Prob. coincidence','FontName','Arial','fontsize',fs)

set(gca,'YTick',yt)
set(gca,'XTick',0:0.1:0.2)
set(gca,'XTickLabel',0:0.1:0.2, 'FontName','Arial','Fontsize',fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,4,3)
hold on
for g=1:ng
    bb=bar(g,msyn(g),0.5,'FaceColor',col{g},'EdgeColor','k','LineWidth',lw,'FaceAlpha',0.5);
end
hold off

box off
grid on
ylim([0,0.075])
xlim([1-.7,ng+.7])

maxy=max(msyn);
ddx=0.1;
pval=pval([1,3,2]);
for ii=1:2
    line([ii+ddx,ii+1-ddx],[maxy+dx,maxy+dx],'color','k')
    line([ii+ddx,ii+ddx],[maxy+dx,maxy+dx/2],'color','k')
    line([ii+1-ddx,ii+1-ddx],[maxy+dx,maxy+ dx/2],'color','k')
    
    if pval(ii)<0.05
        text(ii+.4,maxy+ 1.8*dx,'*','fontsize',fs+3)
    else
        text(ii+.2,maxy+1.8*dx,'n.s.','fontsize',fs)
    end
end

line([1+ddx,3-ddx],[maxy+3*dx,maxy+3*dx],'color','k')
line([1+ddx,1+ddx],[maxy+3*dx,maxy+3*dx-dx/2],'color','k')
line([3-ddx,3-ddx],[maxy+3*dx,maxy+ 3*dx-dx/2],'color','k')

if pval(3)<0.05
    text(1.7,maxy+ 3.5*dx,'*','fontsize',fs+3)
else
    text(1.2,maxy+4*dx,'n.s.','fontsize',fs)
end

title(titles{3}, 'FontName','Arial','Fontsize',fs,'Fontweight','normal')

set(gca,'YTick',[0.02,0.04])
set(gca,'YTickLabel',[0.02,0.04],'FontName','Arial','Fontsize',fs)
set(gca,'XTick',1:ng)
set(gca,'XTickLabel',nameg);

set(gca,'LineWidth',lwa,'TickLength',[0.025 0.025]);

%%%%%%%%%%%%%%%%%%%%%%
subplot(1,4,4)

x=tauvec;
hold on
for ii=1:ng
    y1=mc{ii} - semc{ii};
    y2=mc{ii} + semc{ii};
    patch([x fliplr(x)], [y1' fliplr(y2')], col{ii},'FaceAlpha',0.5,'EdgeColor','k')
end
hold off

axis([0,tauvec(end),0,0.25])
xlabel('maximum lag [ms]','FontName','Arial','Fontsize',fs)
title(titles{4}, 'FontName','Arial','Fontsize',fs,'Fontweight','normal')

set(gca,'XTick',0:10:20)
set(gca,'XTickLabel',0:10:20,'FontName','Arial','Fontsize',fs)
set(gca,'YTick',[0.1,0.2])
set(gca,'YTickLabel',[0.1,0.2],'FontName','Arial','Fontsize',fs)

set(gca,'LineWidth',lwa,'TickLength',[0.025 0.025]);

%%%%%%%%%%%%%


set(H, 'Units','centimeters', 'Position', pos_vec)
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)]) % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%%

