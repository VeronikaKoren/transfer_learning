
clear all
close all
clc

%%

savefig=0;
titles={'all conditions','correct non-match','incorrect non-match','match'};
icase=2;

K=350;

%%
ba=1;
period=2;

namea={'V1','V4'};
namep={'target','test'};
namei={'all','cnm','inm','cm'};
nameg={'SG','G','IG'};
ng=length(nameg);

figname=['cm_layer_inout_',namei{icase},'_',sprintf('%1.0i',K)];
savefile='/home/veronika/Dropbox/transfer/figure/layer/';


green=[0,0.7,0.5];

fs=9;
ms=5;
lw=1.2;
lwa=1;

%% plot cm

addpath('/home/veronika/synced/transfer_result/pairwise/cm/across_pairwise/')
loadname=['cm_layer_across_',namei{icase},'_',sprintf('%1.0i',K),'.mat'];                
load(loadname)

%%
rmat=cell(3,2);
for d=1:2
    rmat{1,d}=cell2mat(cellfun(@(x) x(:,d),a1,'un',0));
    rmat{2,d}=cell2mat(cellfun(@(x) x(:,d),a2,'un',0));
    rmat{3,d}=cell2mat(cellfun(@(x) x(:,d),a3,'un',0));
end

%% define incoming and outgoing links

% outgoing
g_out=cat(1,rmat{1,1}, rmat{3,2});
sg_out=cat(1,rmat{1,2}, rmat{2,2});
ig_out=cat(1,rmat{2,1}, rmat{3,1});
% incoming
g_in=cat(1,rmat{1,2}, rmat{3,1});
sg_in=cat(1,rmat{1,1}, rmat{2,1});
ig_in=cat(1,rmat{2,2}, rmat{3,2});

%% put in a cell and put labels in the same form

inout=cell(6,1);
inout{1}=sg_out;
inout{2}=g_out;
inout{3}=ig_out;
inout{4}=sg_in;
inout{5}=g_in;
inout{6}=ig_in;

nameio=cell(6,1);
nameio{1}='sg_{out}';
nameio{2}='g_{out}';
nameio{3}='ig_{out}';
nameio{4}='sg_{in}';
nameio{5}='g_{in}';
nameio{6}='ig_{in}';

means=cellfun(@mean, inout);
display(means)

%% prepare kruskal-wallis

x_all=[];
name_all=[];

for c=1:6
     
    x=inout{c};
    x_all=cat(1,x_all,x);
    
    replica=repmat(nameio(c),length(x),1);
    name_all=cat(1,name_all,replica);
    
end

%% test kruskal-wallis and multiple comparisons

[p,~,stats]=kruskalwallis(x_all,name_all,'off');
display(p,'pval kruskal-wallis')

z=multcompare(stats,'display','on','ctype','tukey-kramer');
pmc=z(:,[1,2,6]);

%% get p-values

compare={[1,4],[2,5],[3,6]};
pval=zeros(3,1);
for c=1:3
    idxes=find(pmc(:,1)==compare{c}(1));
    subidx=find(pmc(idxes,2)==compare{c}(2));
    index=idxes(subidx);
    pval(c)=pmc(index,3);
end

%% weights

w=cellfun(@mean, inout);
wnorm=5*(w-min(w))/max(w-min(w))+0.5;

[~,idx_outgoing]=max(w(1:3));
[~,idx_incoming]=max(w(4:6));

%%
pos_vec=[0,0,8,13];

H=figure('name',figname,'visible','on');
for ii=1:3
    
    G=digraph(1,2,wnorm(ii+3));
    G = addedge(G,2,3,wnorm(ii));
    wh=G.Edges.Weight;
    
    subplot(3,1,ii)
    
    g=plot(G,'Layout','force','EdgeColor',green);
    if ii==idx_incoming
        rectangle('Position',[g.XData(1) g.YData(1) abs(g.XData(2)-g.XData(1)) abs(g.YData(2)-g.YData(1))],'FaceColor',[1,1,0,0.2],'EdgeColor',[1,0,0,0.5],'Linewidth',1.5)
    elseif ii==idx_outgoing
        rectangle('Position',[g.XData(2) g.YData(2) abs(g.XData(3)-g.XData(2)) abs(g.YData(3)-g.YData(2))],'FaceColor',[1,1,0,0.2],'EdgeColor',[1,0,1],'Linewidth',1.5)
    end
    
    g.LineWidth=wh;
    g.NodeColor='k';
    g.NodeLabel={'in',nameg{ii},'out'};
    g.MarkerSize=6;
    g.LineStyle='-';
    g.ArrowPosition=0.7;
    g.ArrowSize=16;
    g.NodeFontName='Arial';
    g.NodeFontSize=12;
    g.NodeLabelColor=[0,0,0];
    g.NodeFontWeight='normal';
    
    box off
    
    %set(gca,'Visible','off')      
end
axes
h1=title(titles{icase},'FontName','Arial','Fontsize',fs+3,'FontWeight','bold');

set(gca,'Visible','off')
set(h1,'visible','on')
set(H, 'Units','centimeters', 'Position', pos_vec)
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)]) % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

