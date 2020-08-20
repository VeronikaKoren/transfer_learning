% plot cv2

clear all
close all
clc

savefig=1;
ba=1;

addpath('/home/veronika/synced/transfer_result/basic_stat/cv2/')

namep={'target','test'};
namea={'V1','V4'};
namebeh={'different','same'};

figname='cv2';
savefile='/home/veronika/Dropbox/transfer/figure/f1/';


%%

fs=9;
lw=1.5;         % linewidth for plots
lwb=3;          % width for the black bar
lwa=1;

%% load

cv2=cell(2,1);
sem=cell(2,1);
    
for period=1:2
    
    loadname=['cv2_',namea{ba},'_',namep{period}];
    load(loadname)
    cv2{period}=cv_sorted;
    sem{period}=sem_sorted;
    
end

R=cellfun(@(x) corr(x(:,1),x(:,2)),cv2);
%% small plot
%% plot

pos_vec=[0,0,10,6];
ax=[0.5,1.5,0.5,1.5];
H=figure('name',figname);


for p=1:2
    
    subplot(1,2,p)
    scatter(cv2{p}(:,1),cv2{p}(:,2),'.k')
    axis(ax)
    if p==1
        ylabel('cv_2 "same"','FontName','Arial','fontsize',fs)
    end
    set(gca,'XTick',[0.5,1,1.5],'FontName','Arial','fontsize',fs)
    set(gca,'YTick',[1,1.5],'FontName','Arial','fontsize',fs)
    
    title(namep{p},'FontName','Arial','fontsize',fs,'fontweight','normal')
    op=get(gca,'OuterPosition');
    set(gca,'OuterPosition',[op(1) op(2)+0.08 op(3) op(4)-0.08]); % OuterPosition = [left bottom width height]
    
    text(0.1,0.85,['R=' sprintf('%0.4f',R(p))],'units','normalized','FontName','Arial','fontsize',fs)
end


axes

h2 = xlabel ('cv_2 "different"','units','normalized','Position',[0.5,-.04,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h2,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec)                                                      % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])      % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%% plot

%{
blue=[0,0.48,0.74];
green=[0.2,0.6,0]; 
col={blue,green};

pos_vec=[0,0,9,9];
maxy=1.8;

yt=0.5:0.5:1.5;
xt=[80,160];

H=figure('name',figname);

for p=1:2
    
    subplot(2,2,p)
    hold on
    
    x=1:size(cv2{p},1);
    
    for i=1:2
        y1=cv2{p}(:,i) - sem{p}(:,i);
        y2=cv2{p}(:,i) + sem{p}(:,i);
        
        hold on
        patch([x fliplr(x)], [y1' fliplr(y2')], col{i},'FaceAlpha',0.3,'EdgeColor',col{i})
        
    end
    
    hold off
    
    ylim([0,maxy])
    xlim([0,160])
    
    set(gca,'YTick',yt)
    set(gca,'XTick',xt)
    
    set(gca,'XTickLabel',xt, 'FontName','Arial','fontsize',fs)
    if p==1
        set(gca,'YTickLabel',yt, 'FontName','Arial','fontsize',fs)
    else
        set(gca,'YTickLabel',[])
    end
    title(namep{p},'FontName','Arial','fontsize',fs,'fontweight','normal')
    
    if p==1
        ylabel('cv_2','FontName','Arial','fontsize',fs)
    end
    box off
    
    if p==1
        text(0.3,0.78,namebeh{1},'units','normalized','color',col{1},'FontName','Arial','fontsize',fs)
        text(0.3,0.9,namebeh{2},'units','normalized','color',col{2},'FontName','Arial','fontsize',fs)
    end
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
    
end

ax=[0.5,1.5,0.5,1.5];
for p=1:2
    
    subplot(2,2,p+2)
    scatter(cv2{p}(:,1),cv2{p}(:,2),8,'.k')
    axis(ax)
    if p==1
        ylabel('cv_2 "same"','FontName','Arial','fontsize',fs)
    end
    set(gca,'XTick',[0.5,1,1.5],'FontName','Arial','fontsize',fs)
    set(gca,'YTick',[1,1.5],'FontName','Arial','fontsize',fs)
    
    text(0.1,0.85,['R=' sprintf('%0.4f',R(p))],'units','normalized','FontName','Arial','fontsize',fs)
end

axes
h1 = text (0.3,0.48,'neuron index (sorted)','units','normalized','FontName','Arial','fontsize',fs);
h2 = xlabel ('cv_2 "different"','units','normalized','Position',[0.5,-.07,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h1,'visible','on')
set(h2,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec)                                                      % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])      % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%}