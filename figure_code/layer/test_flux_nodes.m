
clear all
close all
clc

%%

savefig=0;
namei={'cnm','inm','cm'};
info_cases=[1,3];

K=350;
alpha=0.05;

%%
ba=1;
period=2;

namea={'V1','V4'};
namep={'target','test'};
namet={'correct v. incorrect','match. non-match','correct match v. incorrect non-match'};
nameg={'SG','G','IG'};
ng=length(nameg);

figname=['flow_compare_',namei{info_cases(1)},'_',namei{info_cases(2)},'_',sprintf('%1.0i',K)];
savefile='/home/veronika/Dropbox/transfer/figure/layer/';

fs=11;
ms=5;
lw=1.2;
lwa=1;

blue=[0,0.48,0.74]; % for correct NM
red=[0.85,0.32,0.1]; % for match
gray=[0.2,0.2,0.2]; % for incorrect non-match

col={red,blue};

%% plot cm

addpath('/home/veronika/transfer_learning/result/pairwise/cm/across/')

flux=cell(2,3);
namef=cell(2,3);
for ic=1:2
    info_case=info_cases(ic);
    
    loadname=['cm_layer_across_',namei{info_case},'_',sprintf('%1.0i',K),'.mat'];
    load(loadname)
    
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
    
    % 
    flux{ic,1}=cat(1,sg_out,sg_in);
    flux{ic,2}=cat(1,g_out,g_in);
    flux{ic,3}=cat(1,ig_out,ig_in);
    
    
    namef{ic,1}=['SG' ,' in ',namei{info_case}];
    namef{ic,2}=['G', ' in ',namei{info_case}];
    namef{ic,3}=['IG',' in ',namei{info_case}];
    
   
end

%% prepare kruskal-wallis

x_all=[];
name_all=[];
for ic=1:2
    for c=1:3
        
        x=flux{ic,c};
        x_all=cat(1,x_all,x);
        
        replica=repmat(namef(ic,c),length(x),1);
        name_all=cat(1,name_all,replica);
        
    end
end
%% test kruskal-wallis and multiple comparisons

[p,~,stats]=kruskalwallis(x_all,name_all,'off');
display(p,'pval kruskal-wallis')

z=multcompare(stats,'display','on','ctype','tukey-kramer');
pmc=z(:,[1,2,6]);

%% get p-values

compare=[[1;4],[2;5],[3;6]];
pval=zeros(3,1);
for c=1:3
    idxes=find(pmc(:,1)==compare(1,c));
    subidx=find(pmc(idxes,2)==compare(2,c));
    index=idxes(subidx);
    pval(c)=pmc(index,3);
end
display(pval,'p-values compare cases')

%% plot flow in layers
if sum(info_cases)==3
    titles=namet{1}
elseif sum(info_cases)==4
    titles=namet{2}
else
    titles=namet{3}
end
xt=0:0.1:0.3;
yt=0:0.5:1;
pos_vec=[0,0,12,13];

%%%%%%%%%%%%%%%%
ylim2=[0.23,0.28];
xb1=1;
xb2=2;
dy=(ylim2(2)-ylim2(1))/20;

dx=(xb2-xb1)/100;
xmax=ylim2(2);
%%%%%%%%%%%%%%%%%%%%%%%

pltidx=[1,3,5];
H=figure('name',figname);

for g=1:3
    
    y1=flux{2,g};
    y2=flux{1,g};
    
    [f1,x1]=ecdf(y1);
    [f2,x2]=ecdf(y2);
     
    subplot(3,2,pltidx(g))
    hold on
    plot(x1,f1,'LineWidth',lw,'color',col{1})
    plot(x2,f2,'LineWidth',lw,'color',col{2})
    
    hold off
    grid on
    yticks(yt)
    xticks(xt)
    
    axis([0.1,0.4,0,1])
    set(gca,'YTickLabel',yt, 'FontName','Arial','fontsize',fs)
    
    idxi=[2,1];
    if g==1
        for ii=1:2
            jj=idxi(ii);
            text(0.02,0.65+(ii-1)*0.2,namei{jj},'units','normalized','FontName','Arial','fontsize',fs,'color',col{ii})
        end
    end
        
    if g==3
        set(gca,'XTickLabel',xt, 'FontName','Arial','fontsize',fs)
        xlabel ('flow','FontName','Arial','fontsize',fs)
    else
        set(gca,'XTickLabel',[], 'FontName','Arial','fontsize',fs)
    end
    
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mflow=cat(2,nanmean(y1),nanmean(y2));
    subplot(3,2,pltidx(g)+1)
    hold on
    for ii=1:2
        bb=bar(ii,mflow(ii),0.5,'FaceColor',col{ii},'EdgeColor','k','LineWidth',lw,'FaceAlpha',0.5);
    end
    hold off
    ylim(ylim2)
    text(0.95,0.5,nameg{g}, 'units','normalized','FontName','Arial','fontsize',fs,'FontWeight','normal') 
    
    set(gca,'XTick',[1,2])
    set(gca,'XTickLabel',[])
    set(gca,'YTick',[0.24,0.26])
    set(gca,'YTickLabel',[0.24,0.26],'FontName','Arial','fontsize',fs)
    
    % significance line  
    xb=max(mflow)+(0.02*max(mflow));
 
    line([xb1,xb2],[xb, xb],'color','k')
    line([xb1,xb1],[xb-dy,xb],'color','k')
    line([xb2,xb2],[xb-dy,xb],'color','k')
    
    if pval(g)<alpha
        th=text((xb1+xb2)/2-0.05*xmax,xb+dy,'*','fontsize',fs+2);
    else
        th=text((xb1+xb2)/2-xmax,xb+2*dy,'n.s.','fontsize',fs-1);
    end
    
    op=get(gca,'OuterPosition');
    set(gca,'OuterPosition',[op(1)+0.02 op(2) op(3)-0.02 op(4)]); % OuterPosition = [left bottom width height]
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
    
end

axes
h0 = title(titles,'FontName','Arial','fontsize',fs);
h2 = ylabel ('empirical Cumulative Distr. Fun.','units','normalized','Position',[-0.1,0.5,0],'FontName','Arial','fontsize',fs);

set(gca,'Visible','off')
set(h2,'visible','on')
set(h0,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec) % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end


%% plot difference of flow between match and non-match in layers
%{
xt=-0.1:0.1:0.1;
yt=0:0.04:0.04;
pos_vec=[0,0,8,13];

ylim=[0,0.08];

x1=-0.025;
x2=0.025;
dy=(ylim(2)-ylim(1))/20;

dx=(x2-x1)/100;
xmax=ylim(2);

H=figure('name',figname);

for g=1:3
    
    [f,xvec]=ksdensity(flux{2,g}-flux{1,g});
    m=nanmean(flux{2,g}-flux{1,g});
    fnorm=f./sum(f);
    
    subplot(3,1,g)
    hold on
    plot(xvec,fnorm,'LineWidth',lw,'color','k')
    hm=plot(m,max(fnorm),'d','markersize',ms,'MarkerFaceColor','r');
    hold off
    
    line([0,0],[0,0.055],'color',[0.5,0.5,0.5])
    axis([-0.15,0.15,ylim(1),ylim(2)])
    
    grid on
    yticks(yt)
    xticks(xt)
    
    set(gca,'XTickLabel',xt, 'FontName','Arial','fontsize',fs)
    set(gca,'YTickLabel',yt, 'FontName','Arial','fontsize',fs)
    
    text(0.95,0.5,nameg{g}, 'units','normalized','FontName','Arial','fontsize',fs,'FontWeight','normal')
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
    
    % significance line
       
    xb=max(fnorm)+(0.25*max(fnorm));
 
    line([x1,x2],[xb, xb],'color','k')
    line([x1,x1],[xb-dy,xb],'color','k')
    line([x2,x2],[xb-dy,xb],'color','k')
    
    if pval (g)<alpha
        th=text((x1+x2)/2-0.05*xmax,xb+dy,'*','fontsize',fs+2);
    else
        th=text((x1+x2)/2-0.12*xmax,xb+2*dy,'n.s.','fontsize',fs-1);
    end
    
    op=get(gca,'OuterPosition');
    set(gca,'OuterPosition',[op(1)+0.08 op(2) op(3)-0.08 op(4)]); % OuterPosition = [left bottom width height]
    
end

axes
h1 = xlabel (['total flow ',titles{info_cases(2)},' - ', titles{info_cases(2)}],'units','normalized','Position',[0.5,-0.07,0],'FontName','Arial','fontsize',fs);
h2 = ylabel ('probability density','units','normalized','Position',[-0.08,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h2,'visible','on')
set(h1,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec) % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end
%}