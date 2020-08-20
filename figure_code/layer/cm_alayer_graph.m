
clear all
close all
clc

%%

savefig=0;

titles={'correct non-match','incorrect non-match',' match'};
info_case=1;
display(['plot directed graph ',titles{info_case}])

K=350;

%%  

ba=1;
period=2;

namea={'V1','V4'};
namep={'target','test'};
namei={'cnm','inm','cm'};
nameg={'SG','G','IG'};
ng=length(nameg);

figname=['layer_nodes_',namei{info_case},'_',sprintf('%1.0i',K)];
savefile='/home/veronika/Dropbox/transfer/figure/layer/';

green=[0,0.7,0.5];
pos_vec=[0,0,10,8];

fs=9;
ms=5;
lw=1.2;
lwa=1;

%% add path and load result

addpath('/home/veronika/synced/transfer_result/pairwise/cm/across/')
loadname=['cm_layer_across_',namei{info_case},'_',sprintf('%1.0i',K),'.mat'];                
load(loadname)

%% put results in a matrix

rmat=cell(3,2);
for d=1:2
    rmat{1,d}=cell2mat(cellfun(@(x) x(:,d),a1,'un',0));
    rmat{2,d}=cell2mat(cellfun(@(x) x(:,d),a2,'un',0));
    rmat{3,d}=cell2mat(cellfun(@(x) x(:,d),a3,'un',0));
end

%% kruskal-wallis

x_all=[];
name_all=[];
namec=[];
r6pack=[];
for c=1:3
    for dir=1:2
        
        x=rmat{c,dir};
        x_all=cat(1,x_all,x);
        
        replica=repmat(labels(c,dir),length(x),1);
        name_all=cat(1,name_all,replica);
        
        namec=cat(1,namec,labels(c,dir));
        r6pack=cat(1,r6pack,rmat(c,dir));
    end
end

[p,~,stats]=kruskalwallis(x_all,name_all,'off');
display(p,'pval kruskal-wallis')

z=multcompare(stats,'display','on','ctype','tukey-kramer');
pmc=z(:,[1,2,6]);

%% get p-values

compare={[2,4],[1,6],[3,5]};
pval=zeros(3,1);
for c=1:3
    idxes=find(pmc(:,1)==compare{c}(1));
    subidx=find(pmc(idxes,2)==compare{c}(2));
    index=idxes(subidx);
    pval(c)=pmc(index,3);
end

%% ecdf 

figure('units','centimeters','Position',[0,0,12,21],'visible','off')

for c=1:3
    subplot(3,1,c)
    hold on
    ecdf(r6pack{compare{c}(1)})
    ecdf(r6pack{compare{c}(2)})
    
    hold off
    text(0.2,0.8,namec{compare{c}(1)},'units','normalized','color','b')
    text(0.2,0.7,namec{compare{c}(2)},'units','normalized','color','r')
    text(0.15,0.5,['p = ' sprintf('%0.3f',pval(c))],'units','normalized','color','k')
   
    if pval(c)<0.05
        text(0.1,0.47,'*','units','normalized','color','k','fontsize',15)
    end
    xlim([0,0.6])
    
end

%% build adjaceny matrix

means=cellfun(@(x) nanmean(x), rmat);
A=zeros(3,3);
A(1,2)=means(1,2); % SG to G
A(1,3)=means(2,2); % SG to IG
A(2,1)=means(1,1); % G to SG
A(2,3)= means(3,2); % G to IG
A(3,1)=means(2,1);  % IG to SG
A(3,2)= means(3,1); % IG to G


display(labels,'labels')
display(means)

%% make graph with links between layers

G=digraph(A);
w=G.Edges.Weight;
wnorm=5*(w-min(w))/max(w-min(w))+0.5; % normalization

%% plot graph

H=figure('name',figname);
g=plot(G,'r','EdgeColor',green,'LineWidth',wnorm,'NodeLabel',{'SG','G','IG'});

g.MarkerSize=7;
g.LineStyle='-.';
g.ArrowSize=20;
g.NodeFontName='Arial';
g.NodeFontSize=13;
g.NodeLabelColor=[0,0,0];
g.NodeFontWeight='bold';
set(gca,'Color',[1,1,0,0.15])
box on
title(titles{info_case},'FontName','Arial','Fontsize',fs+3,'FontWeight','bold')

set(H, 'Units','centimeters', 'Position', pos_vec) % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

