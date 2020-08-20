% plots LDR of models with permuted spike timing, whole window or before/after the switch 

clear all 
close all
clc 

savefig=0;

perms={'entire_window','half_window'};

alpha=0.025;       % 2.5 percent on each side

%%

savefile='/home/veronika/Dropbox/transfer/figure/permutations/';
figname='perm_timing';

namebeh={'different','same'};

fs=10; 
lw=1.5;

gray=[0.2,0.2,0.2];
blue=[0,0.48,0.74];
green=[0.2,0.6,0];
col={blue,green};

factor=100;
pos_vec=[0,0,9,12];
%% load the permuted       

addpath('/home/veronika/transfer_learning/result/signal/permute_timing/')

np=length(perms);

xsame=cell(np,1);
xdiff=cell(np,1);

for p=1:np
    
    loadname=[perms{p}];
    load(loadname)
    
    xsame{p}=squeeze(nanmean(xp_same)).*factor;
    xdiff{p}=squeeze(nanmean(xp_diff)).*factor;
    clear xp_same xp_diff
    
end

nperm=size(xsame{1},1);

%% compute the boundaries (leaving alpha percent on the top and on the bottom)

K=size(xsame{1},2);
idx=max(1,floor(alpha/K*nperm));

lower_bound=cell(np,1);
upper_bound=cell(np,1);

for p=1:np
    
    lb=zeros(2,K);
    ub=zeros(2,K);
    for k=1:K
        
        x=sort(xsame{p}(:,k));
        lb(1,k)=x(idx+1);
        ub(1,k)=x(end-idx);
        
        y=sort(xdiff{p}(:,k));
        lb(2,k)=y(idx+1);
        ub(2,k)=y(end-idx);
        
    end
    lower_bound{p}=lb;
    upper_bound{p}=ub;
    
end

%% plot

idx_switch=139;
xarrow=[idx_switch+1,idx_switch+1];                                                 % switch of the population signal
yarrow=-[0.2,0.4];

titles={'whole window','before/after zero crossing'};

maxy=0.75;
yt=[-0.3,0,0.3];
xt=0:100:400;
x=1:K;

H=figure('name',figname);
for p=1:np
    
    y1=squeeze(lower_bound{p}(1,:));
    y2=squeeze(upper_bound{p}(1,:));
    
    z1=squeeze(lower_bound{p}(2,:));
    z2=squeeze(upper_bound{p}(2,:));
    
    subplot(2,1,p)
    hold on
    
    patch([x fliplr(x)], [y1 fliplr(y2)], col{2},'FaceAlpha',0.3,'EdgeColor', col{2})
    patch([x fliplr(x)], [z1 fliplr(z2)], col{1},'FaceAlpha',0.3,'EdgeColor', col{1})
    plot(1:K,zeros(K,1),'--','color',gray,'linewidth',lw)
    
    if p==2
        line(xarrow-0.5,yarrow,'linewidth',lw+1,'color','r')
        plot(xarrow(1),yarrow(1),'^','markersize',8,'MarkerFaceColor','red','markeredgecolor','red')
        %annotation('arrow',xarrow,yarrow,'color','red')
        text(xarrow(1)-30,yarrow(1)-0.28,'crossing','units','data','color','red','fontsize',fs,'fontname','Arial')
    end
    
    hold off
    box off
    
    ylim([-maxy+0.1,maxy])
    xlim([-2,K+2])
    
    set(gca,'XTick',xt)
    set(gca,'XTickLabel',[])
    
    set(gca,'XTickLabel',xt, 'FontName','Arial','fontsize',fs)
 
    set(gca,'YTick',yt)
    set(gca,'YTickLabel',yt, 'FontName','Arial','fontsize',fs)
    
    %{
    if p==1
        text(0.1,0.87, 'different', 'color', col{1},'units','normalized', 'FontName','Arial','fontsize',fs)
        text(0.1,0.95, 'same', 'color', col{2},'units','normalized', 'FontName','Arial','fontsize',fs)
    end
    %}
    
    title(titles{p},'Fontweight','normal','FontName','Arial','fontsize',fs)
    set(gca,'LineWidth',1.2,'TickLength',[0.025 0.025]);
    
    op=get(gca,'OuterPosition');
    set(gca,'OuterPosition',[op(1)-0.02 op(2)+0.03 op(3)+0.02 op(4)-0.03]); % OuterPosition = [left bottom width height]
    %xlabel ('time (ms)','FontName','Arial','fontsize',fs)
end


axes

h1 = xlabel ('time (ms)','units','normalized','Position',[0.5,-0.06,0],'FontName','Arial','fontsize',fs);
%h2 = ylabel ('population signal (a.u.)','units','normalized','Position',[-0.11,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
%set(h2,'visible','on')
set(h1,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec) % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)]) % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

