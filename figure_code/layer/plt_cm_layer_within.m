
clear all
close all
clc

savefig=0;
namei={'s+c','c','all','cnm','inm','cm'};
info_case=3;

K=500;

%%
ba=1;
period=2;

namea={'V1','V4'};
namep={'target','test'};
nameg={'SG','G','IG'};
ng=length(nameg);

figname=['cm_layer_within_',namei{info_case},'_',sprintf('%1.0i',K)];
savefile='/home/veronika/Dropbox/transfer/figure/layer/';

pos_vec=[0,0,12,5];

fs=10;
ms=5;
lw=1.2;
lwa=1;

blue=[0,0.48,0.74];
magenta=[1,0,1];
black=[0,0,0];
col={magenta,blue, black};

%% load

addpath('/home/veronika/synced/transfer_result/pairwise/cm/within/')
loadname=['cm_layer_within_',namei{info_case},'_',sprintf('%1.0i',K),'.mat'];                
load(loadname)

links=cellfun(@(x) cat(1,x(:,1),x(:,2)), cm_lay,'un',0);    % get all links (both directions)

%% kruskal-wallis

x_all=[];
name_all=[];

for g=1:3
    
    x=links{g};
    x_all=cat(1,x_all,x);
    
    replica=repmat(nameg(g),length(x),1);
    name_all=cat(1,name_all,replica);
end

[p,~,stats]=kruskalwallis(x_all,name_all,'off');
display(p,'pval kruskal-wallis')
z=multcompare(stats,'display','off','ctype','tukey-kramer');
pmc=z(:,[1,2,6]);
display(pmc,'multiple comparison test')

%%
ml=cellfun(@mean,links);
seml=cellfun(@(x) std(x)./sqrt(length(x)),links);
pval=pmc([1,3,2],3);


%% plot ecdf and bars

maxy=max(ml);
ylimit=[0.22,0.38];
dy=(ylimit(2)-ylimit(1))/10;
a=0.55;
ddx=0.06;

H=figure('name',figname,'visible','on');

subplot(1,2,1)
yt=0:0.5:1;
xt=0:0.2:0.6;

hold on
for g=1:3
    [f,xi]=ecdf(links{g});
    plot(xi,f,'color',col{g},'linewidth',lw);
    text(0.15,0.75-(g-1)*0.1,nameg{g},'units','normalized','FontName','Arial','fontsize',fs,'color',col{g})
end

xlabel('cross-mapping coeff.', 'FontName','Arial','Fontsize',fs)
ylabel('empirical CDF', 'FontName','Arial','Fontsize',fs)
axis([0,0.6,0,1.1])
grid on

set(gca,'YTick',yt)
set(gca,'YTickLabel',yt, 'FontName','Arial','Fontsize',fs)
set(gca,'XTick',xt)
set(gca,'XTickLabel',xt,'FontName','Arial','Fontsize',fs)

set(gca,'LineWidth',lwa,'TickLength',[0.025 0.025]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
hold on
for g=1:3
    bar(g,ml(g),0.5,'FaceColor',col{g},'EdgeColor','k','FaceAlpha',0.5,'linewidth',lw);
end
hold off

ylabel('cross-mapping', 'FontName','Arial','Fontsize',fs)
ylim(ylimit)
xlim([1-.7,3+.7])
set(gca,'XTick',1:3)
set(gca,'XTickLabel',nameg,'FontName','Arial','Fontsize',fs)


for ii=1:2
    line([ii+ddx,ii+1-ddx],[maxy+dy,maxy+dy],'color','k')
    line([ii+ddx,ii+ddx],[maxy+dy,maxy+ (a*dy)],'color','k')
    line([ii+1-ddx,ii+1-ddx],[maxy+dy,maxy+ (a*dy)],'color','k')
    
    if pval(ii)<0.05
        text(ii+.4,maxy+ 1.3*dy,'*','fontsize',fs+3)
    else
        text(ii+.2,maxy+1.5*dy,'n.s.','fontsize',fs)
    end
end


line([1+2*ddx,3-2*ddx],[maxy+2.5*dy,maxy+2.5*dy],'color','k')
line([1+2*ddx,1+2*ddx],[maxy+2.5*dy,maxy+2.5*dy-(a*dy)],'color','k')
line([3-2*ddx,3-2*ddx],[maxy+2.5*dy,maxy+ 2.5*dy-(a*dy)],'color','k')
if pval(3)<0.05
    text(1.9,maxy+ 3*dy,'*','fontsize',fs+3)
else
    text(1.6,maxy+3*dy,'n.s.','fontsize',fs)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(H, 'Units','centimeters', 'Position', pos_vec)
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)]) % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%%
pval=pmc(:,3);
display(pval,'p-value multcompare layers')

%{
write p-values in the bar plot :)
xtxt=(pmc(:,1)+pmc(:,2))/2-0.3;
ytxt=[maxy,maxy+0.03,maxy];
for c=1:3
    text(xtxt(c),ytxt(c),['p = ' sprintf('%0.3f',pval(c))])
end
%}