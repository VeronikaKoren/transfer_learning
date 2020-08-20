% plot average firing rates

clear all
close all
clc

savefig=0;
ba=1;

addpath('/home/veronika/synced/transfer_result/basic_stat/frate/')
addpath('/home/veronika/synced/transfer_result/basic_stat/cv2/')

namep={'target','test'};
namea={'V1','V4'};
namebeh={'different','same'};

figname='frate_cv';
savefile='/home/veronika/Dropbox/transfer/figure/f1/';
pos_vec=[0,0,9,9];

%%
fs=9;
lw=1.5;         % linewidth for plots
lwb=3;          % width for the black bar
lwa=1;

blue=[0,0.48,0.74];
green=[0.2,0.6,0]; 
col={blue,green};

%% load f.rate

frate=cell(2,1);
sem_rate=cell(2,1);
    
for period=1:2
    
    loadname=['frate_',namea{ba},'_',namep{period}];
    load(loadname)
    frate{period}=fr_sorted;
    sem_rate{period}=sem_sorted;
    
end

R_rate=cellfun(@(x) corr(x(:,1),x(:,2)),frate);

%% load cv2

cv2=cell(2,1);
sem_cv=cell(2,1);
    
for period=1:2
    
    loadname=['cv2_',namea{ba},'_',namep{period}];
    load(loadname)
    cv2{period}=cv_sorted;
    sem_cv{period}=sem_sorted;
    
end

R_cv=cellfun(@(x) corr(x(:,1),x(:,2)),cv2);
%% plot

maxy=80;

yt=0:30:60;
xt=[80,160];

H=figure('name',figname);

for p=1:2
    subplot(2,2,p)
    hold on
    scatter(frate{p}(:,1),frate{p}(:,2),'.k')
    if p==1
        ylabel('spikes/sec "same"','FontName','Arial','fontsize',fs)
    end
    set(gca,'XTick',[0,50,100],'FontName','Arial','fontsize',fs)
    set(gca,'YTick',[0,50,100],'FontName','Arial','fontsize',fs)
    
    text(0.1,0.85,['R=' sprintf('%0.4f',R_rate(p))],'units','normalized','FontName','Arial','fontsize',fs)
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
    
    text(0.1,0.85,['R=' sprintf('%0.4f',R_cv(p))],'units','normalized','FontName','Arial','fontsize',fs)
    
end

axes
h1 = xlabel ('cv_2 "different"','units','normalized','Position',[0.5,-.06,0],'FontName','Arial','fontsize',fs);
h2 = xlabel ('spikes/sec "different"','units','normalized','Position',[0.5,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h2,'visible','on')
set(h1,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec)                                                      % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])      % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%%

mean_tar=mean(frate{1}(:))
mean_test=mean(frate{2}(:))

