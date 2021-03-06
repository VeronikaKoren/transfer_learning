
clear all
close all
clc

savefig=0;

ba=1;
period=2;
K=400;

%%
col={'m','k'};

lw=1.0;     % linewidth
ms=4;       % markersize
fs=10;      % fontsize
lwa=1;

namea={'V1','V4'};
namep={'target','test'};
nameg={'bursty','non-bursty'};

figname=['ccg_bursty_',sprintf('%1.0i',K)];
savefile='/home/veronika/Dropbox/transfer/figure/bursty/';

%% plot ccg sign for plus and minus

addpath('/home/veronika/transfer_learning/result/pairwise/ccg/')

loadname=['ccg_bursty_',sprintf('%1.0i',K),'.mat'];                
load(loadname)
    
%%
variable=cell(2,1);
variable{1}=cell2mat(g1);
variable{2}=cell2mat(g2);

%%

iW=(size(variable{1,1},2)+1)/2;
speak=cellfun(@(x) x(:,iW),variable,'UniformOutput',false);
display(pval,'p-value permutation test');

xvec=-iW+1:iW-1;
msyn=cellfun(@mean,speak);
semsyn=cellfun(@(x) std(x)./sqrt(length(x)),speak);

%%
coeffs=cell(2,1);
coeffs{1}=r1'; coeffs{2}=r2';
mc=cellfun(@mean, coeffs,'UniformOutput', false);
semc=cellfun(@(x) std(x)./sqrt(size(x,1)),coeffs,'UniformOutput', false);

%% plot

pos_vec=[0 0 18 5];

titles={'Cross-corr.','Synchrony','Average synchrony','Area under ccg'};

yt=0:0.04:0.12;
ylimit=[-0.02,0.1];
dx=abs(ylimit(2)-ylimit(1))/20;

H=figure('name',figname,'visible','on');
%%%%%%%%%%%%%%%%%%%

[~,idx]=sort(speak{1});

idx1=idx(2*round(length(idx)/3));
[~,idx]=sort(speak{2});
idx2=idx(2*round(length(idx)/3));

subplot(1,4,1)
xt=[-100,0,100];
hold on
plot(xvec,variable{1}(idx1,:),'color',[col{1},0.5],'linewidth',lw+0.5);
plot(xvec,variable{2}(idx2,:),'color',[col{2},0.5],'linewidth',lw);
hold off

for ii=1:2
    text(0.5,0.9-(ii-1)*0.1,nameg{ii},'units','normalized','FontName','Arial','fontsize',fs,'color',col{ii})
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yt=[0,0.5,1];
[f,x1]=ecdf(speak{1});
[g,x2]=ecdf(speak{2});

subplot(1,4,2)
hold on
for ii=1:2
    [f,xi]=ecdf(speak{ii});
    plot(xi,f,'color',col{ii},'linewidth',lw)
end

hold off

grid on
xlim([-0.01,0.2])

title(titles{2}, 'FontName','Arial','Fontsize',fs,'Fontweight','normal')
xlabel ('Prob. coincidence','FontName','Arial','fontsize',fs)

set(gca,'YTick',yt)
set(gca,'XTick',0:0.1:0.2)
set(gca,'XTickLabel',0:0.1:0.2, 'FontName','Arial','Fontsize',fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,4,3)
hold on
for b=1:2
    bb=bar(b,msyn(b),0.5,'FaceColor',col{b},'EdgeColor','k','LineWidth',lw,'FaceAlpha',0.5);
end
hold off

box off
grid on
ylim([0,0.06])
xlim([0.3,2.7])

maxy=max(msyn);

line([1,2],[maxy+dx,maxy+dx],'color','k')
line([1,1],[maxy+dx,maxy+dx/2],'color','k')
line([2,2],[maxy+dx,maxy+ dx/2],'color','k')

if pval<0.05
    text(1.4,maxy+ 1.8*dx,'*','fontsize',fs+3)
else
    text(1.1,maxy+1.8*dx,'n.s.','fontsize',fs)
end

title(titles{3}, 'FontName','Arial','Fontsize',fs,'Fontweight','normal')

set(gca,'YTick',[0.02,0.04])
set(gca,'YTickLabel',[0.02,0.04],'FontName','Arial','Fontsize',fs)
set(gca,'XTick',[1,2])
set(gca,'XTickLabel',[]);

set(gca,'LineWidth',lwa,'TickLength',[0.025 0.025]);

%%%%%%%%%%%%%%%%%%%%%%
subplot(1,4,4)

x=tauvec;
hold on
for ii=1:2
    y1=mc{ii} - semc{ii};
    y2=mc{ii} + semc{ii};
    patch([x fliplr(x)], [y1 fliplr(y2)], col{ii},'FaceAlpha',0.5,'EdgeColor','k')
end
hold off

axis([0,tauvec(end),0,0.18])
xlabel('maximum lag [ms]','FontName','Arial','Fontsize',fs)
ylabel('corr. coeffcient','FontName','Arial','Fontsize',fs)
title(titles{4}, 'FontName','Arial','Fontsize',fs,'Fontweight','normal')

set(gca,'XTick',0:10:20)
set(gca,'XTickLabel',0:10:20,'FontName','Arial','Fontsize',fs)
set(gca,'YTick',[0.06,0.12])
set(gca,'YTickLabel',[0.06,0.12],'FontName','Arial','Fontsize',fs)

set(gca,'LineWidth',lwa,'TickLength',[0.025 0.025]);

%%%%%%%%%%%%%

set(H, 'Units','centimeters', 'Position', pos_vec)
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)]) % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
