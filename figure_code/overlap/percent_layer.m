
clear all
close all
clc
                                                                  
K=400;
savefig=0;
info_case=1;

%%

figname='percent_layers';
savefile='/home/veronika/Dropbox/transfer/figure/weights/';

ba=1;     
period=2;

namep={'tar','test'};
namea={'V1','V4'};
namelay={'SG','G','IG'};
namei={'s+c','c'};

%%

addpath('/home/veronika/transfer_learning/result/weight/tags/')
addpath('/home/veronika/transfer_learning/result/weight/weight_alltr/')

loadname=['tag_info_',namei{info_case},'_',sprintf('%1.0i',K)];
load(loadname)

loadname1=['weight_alltr_',namei{info_case},'_',namea{ba},'_', sprintf('%1.0i',K)];
load(loadname1)
wplus=cellfun(@(x) sign(x)==1,w_alltr,'un',0);

loadname3=['ncell_lay_',namea{ba}];
load(loadname3);

loadname4=['tag_bursty_',sprintf('%1.0i',K)];
load(loadname4)

nbses=length(tag_info);
%%
ninfo=zeros(3,1);
for r=1:3
    
    nb_info=0;
    for sess=1:nbses
        
        ncl=ncell_lay(sess,:);
        s=cumsum([0,ncl]);
        delta = s(r) + 1 : s(r+1);
        
        nb=sum(tag_info{sess}(delta));    % count the number of plus neurons in each layer
        nb_info=nb_info+nb;
        
    end
    
    ninfo(r)=nb_info;
    
end

percent_info=ninfo'./sum(ncell_lay);     % normalize with the number of neurons per layer
display(percent_info,'percent info in layers')
%%

nsign=zeros(3,1);
for r=1:3
    
    nb_sign=0;
    for sess=1:nbses
        
        ncl=ncell_lay(sess,:);
        s=cumsum([0,ncl]);
        delta = s(r) + 1 : s(r+1);
        
        nb=sum(wplus{sess}(delta));
        nb_sign=nb_sign+nb;
        
    end
    
    nsign(r)=nb_sign;
    
end

percent_sign_lay=nsign'./sum(ncell_lay);
display(percent_sign_lay,'percent sign in layers')


%%

nb_bursty=zeros(3,1);
for r=1:3
    
    nb_bur=0;
    for sess=1:nbses
        
        ncl=ncell_lay(sess,:);
        s=cumsum([0,ncl]);
        delta = s(r) + 1 : s(r+1);
        
        nb=sum(bursty{sess}(delta));
        nb_bur=nb_bur+nb;
        
    end
    
    nb_bursty(r)=nb_bur;
    
end

percent_bur_lay=nb_bursty'./sum(ncell_lay);
display(percent_bur_lay,'percent bursty in layers')

%%

fs=10; % figure settings
lw=1.5;
pos_vec=[0,0,9,6];

pltx=cell(2,1);
pltx{1}=percent_sign_lay;
pltx{2}=percent_bur_lay;

%%

titles={'plus neurons','bursty neurons'};
col=[0,0,0]; % black
yt=0.4:0.05:0.6;

H=figure('name',figname,'visible','on');
for ii=1:2
    subplot(1,2,ii)
    
    bar(pltx{ii},'FaceColor',col,'FaceAlpha',0.3,'EdgeColor',col)
    grid on
    line([0,4],[0.5,0.5],'color','k','linestyle','--')
    ylim([0.35,0.65])
    
    set(gca,'YTick',yt)
    if ii==1
       set(gca,'YTickLabel',yt, 'FontName','Arial','fontsize',fs)
       
    else
        set(gca,'YTickLabel',[])
    end
    set(gca,'XTickLabel',namelay,'FontName','Arial','fontsize',fs)
    title(titles{ii},'FontName','Arial','fontsize',fs,'fontweight','normal')
    box off
    
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
    
    op=get(gca,'OuterPosition');
    set(gca,'OuterPosition',[op(1) op(2)+0.08 op(3) op(4)-0.08]); % OuterPosition = [left bottom width height]
    
end

axes
%text(-0.15,1.05,'A','units','normalized','FontName','Arial','fontsize',fs,'fontweight','bold')
h1 = ylabel ('proportion','units','normalized','Position',[-0.11,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h1,'visible','on')

set(H, 'Units','centimeters', 'Position', pos_vec)                                                      % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])      % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end


