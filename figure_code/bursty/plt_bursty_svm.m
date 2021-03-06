% plot balanced accuracy in sessions and average across sessions

format long

clear all
close all
clc

savefig=0;
%%

T=400;

ba=1;
period=2;

namep={'target','test'};
namea={'V1','V4'};
nameg={'bursty','non-bursty'};

savefile='/home/veronika/Dropbox/transfer/figure/bursty/';
figname=['bac_bursty_',sprintf('%1.0i',T)];

lw=1.0;     % linewidth
ms=4;       % markersize
fs=10;      % fontsize
lwa=1;

blue=[0,0.48,0.74];
red=[0.85,0.32,0.1];
col={'m','k'};

%% load the BAC of the linear SVM for sign

addpath('/home/veronika/transfer_learning/result/weight/bac/')

loadname=['bac_bursty_',sprintf('%1.0i',T),'.mat'];
load(loadname);
    
nbses=length(bac_group);
nperm=size(bac_group{1},2);

% average across permutations
mbac=cell2mat(cellfun(@(x) nanmean(x,2), bac_group, 'UniformOutput', false)');

%% test with t-test BAC for groups are significantly different

bacg=zeros(2,nperm);
for ii=1:2
    bacg(ii,:)=nanmean(cell2mat(cellfun(@(x) x(ii,:),bac_group,'UniformOutput', false)));
end
%%

[h,pval]=ttest2(bacg(1,:),bacg(2,:)); % average across sessions, then test across permutations
display(pval,'p-value difference accuracy bursty/nonbursty')

%% plot sessions and average

[val,sorder]=sort(mbac(1,:)); % session order
order=flip(sorder);

pos_vec=[0,0,15,6];
yt=0.5:0.1:0.7; 
xt=5:5:20;
yt2=[0.5,0.55];

%%%%%%%%

H=figure('name',figname,'visible','on');

subplot(1,3,[1,2])
hold on

plot(0:nbses+1,ones(nbses+2,1).*0.5,'--','color',[0.5,0.5,0.5,0.5],'linewidth',lw)
plot(1:nbses,mbac(1,order),'x','color',col{1},'markersize',ms,'Linewidth',lw);
plot(1:nbses,mbac(2,order),'+','color',col{2},'markersize',ms,'Linewidth',lw);
hold off

text(0.5,0.8,'bursty','units','normalized','color',col{1},'fontsize',fs,'FontName','Arial')
text(0.5,0.7,'non-bursty','units','normalized','color',col{2},'fontsize',fs,'FontName','Arial')
   
axis([0,nbses+1,0.45,0.75])
box off
grid on

set(gca,'XTick',xt)
set(gca,'XTickLabel',xt,'FontName','Arial','fontsize',fs)

set(gca,'YTick',yt)
set(gca,'YTickLabel',yt,'fontsize',fs,'FontName','Arial')
set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
xlabel ('session index (sorted)','FontName','Arial','fontsize',fs);
ylabel ('balanced accuracy','FontName','Arial','fontsize',fs);

op=get(gca,'OuterPosition');
set(gca,'OuterPosition',[op(1) op(2)+0.05 op(3) op(4)-0.05]); % OuterPosition = [left bottom width height]

%%%%%%%%%

bars=mean(mbac,2);

subplot(1,3,3)
hold on
for b=1:2
    bb=bar(b,bars(b),0.5,'FaceColor',col{b},'EdgeColor','k','LineWidth',1,'FaceAlpha',0.5);
end
hold off

% significance line
alpha=0.05;
ylim2=[0.48,0.58];
x1=1; 
x2=2; 
xb=max(bars)+0.25*(max(bars)-0.5);
dy=(ylim2(2)-ylim2(1))/20;
dx=(x2-x1)/100;

xmax=2.7;
line([x1,x2],[xb, xb],'color','k')
line([x1,x1],[xb-dy,xb],'color','k')
line([x2,x2],[xb-dy,xb],'color','k')

if pval<alpha
    th=text((x1+x2)/2-0.05*xmax,xb+dy,'*','fontsize',fs+3);
else
    th=text((x1+x2)/2-0.12*xmax,xb+dy,'n.s.','fontsize',fs);
end
%}
xlim([0.3,2.7])
ylim(ylim2)

box off
grid on

set(gca,'XTick',[1,2])
set(gca,'YTick',yt2)
set(gca,'XTickLabel', nameg, 'FontName','Arial','fontsize',fs)
xtickangle(25)
set(gca,'YTickLabel',yt2,'fontsize',fs,'FontName','Arial')

set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);

set(H, 'Units','centimeters', 'Position', pos_vec) 
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)]) 

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end



