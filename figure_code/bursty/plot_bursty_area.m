
close all
clear all
clc

figname='area_ps';
figname2='power_spectrum';

savefig=0;
savefig2=0;

K=400;          % length of the time window
info_case=1;
%%

ba=1;
period=2;

namep={'tar','test'};
namea={'V1','V4'};
namei={'s+c','c'};                                                                  

savefile='/home/veronika/Dropbox/transfer/figure/bursty/';
fs=10; 
lw=1.5;

%% load data 

addpath('/home/veronika/synced/transfer_result/basic_stat/ps/')
loadname=['ps_',namea{ba},namep{period},'_',sprintf('%1.0i',K)];
load(loadname)

loadname2=['ps_poisson_',namea{ba},namep{period},'_',sprintf('%1.0i',K)];
load(loadname2)

%%
addpath('/home/veronika/synced/transfer_result/weight/tags/')
loadname3=['tag_bursty_',sprintf('%1.0i',K)];
load(loadname3)

nf=length(f);   % number of frequencies evauated
nbses=size(psp,1);
nperm=size(psp{1},1);

%% plot area under PS regular and poisson

mata=cell2mat(aps);
low=cell2mat(lower_bound);
tag=cell2mat(bursty);

[val,order]=sort(mata);  % sort area under PS from smallest to biggest
Ntot=length(mata);

%% plot area under PS

col={'k',[0,0,1], 'r'};
pos_vec=[0,0,8,6];

xt=0:50:150;
yt=0.8:0.1:1;

H=figure('name',figname,'visible','on');
hold on
plot(mata(order),'x','color',col{1})
plot(low(order),'-','color',col{2})
plot(tag(order)*0.75,'*','markersize',2,'color',col{3})
hold off

text(0.1,0.9,' area under PS','units','normalized','color','k','FontName','Arial','fontsize',fs)
text(0.1,0.83,'lower bound permuted','units','normalized','color',col{2},'FontName','Arial','fontsize',fs)
%text(0.1,0.76,'significant neuron','units','normalized','color',col{3},'FontName','Arial','fontsize',fs)

grid on
ylim([0.7,1.15])
xlim([0,Ntot])
xlabel('neuron index (sorted)')
ylabel('area under p.spectrum')

set(gca,'XTick',xt)
set(gca,'YTick',yt)
set(gca,'XTickLabel',xt)
set(gca,'YTickLabel',yt)

set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);

set(H, 'Units','centimeters', 'Position', pos_vec)                                                      % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec(3), pos_vec(4)])      % for saving in the right size

if savefig==1
    saveas(H,[savefile,figname],'pdf');
end

%% compute lower bound for renewal power spectum (as a function of frequency)

alpha=0.05;
idx=nperm*alpha;

upper=cell(nbses,1);
lower=cell(nbses,1);

for sess=1:nbses
    
    N=size(psp{sess},2);
    ub=zeros(N,nf);
    lb=zeros(N,nf);
    for ii=1:N
        
        y=squeeze(psp{sess}(:,ii,:));
        
        for fr=1:size(y,2)
            sorted=sort(y(:,fr)); % sort from small to big;
            lb(ii,fr)=sorted(idx);
            ub(ii,fr)=sorted(end-idx);
        end
        
    end
    
    upper{sess}=ub;
    lower{sess}=lb;
    
end


%% plot ps 4 neurons
fs=fs-2.5;
cs=5;                                  % choose a session
vec=5:nf;
x=f(vec);
pos_vec2=[0,0,14,5];

H=figure('name',figname2,'visible','on');
for ii=1:4
    
    z1=lower{cs}(ii,vec);
    z2=upper{cs}(ii,vec);
   
    subplot(1,4,ii)
    hold on
    
    patch([x fliplr(x)], [z1 fliplr(z2)], col{2},'FaceAlpha',0.3,'EdgeColor',col{2},'EdgeAlpha',0.3,'Linewidth',0.5)
    plot(x,ps_sess{cs}(ii,vec),col{1},'linewidth',lw-0.5)
    hold off
    
    xlim([0,500])
    ylim([0.6,1.4])
    
    grid on
    
    set(gca,'XTick',0:250:500)
    set(gca,'XTickLabel',0:250:500, 'FontName','Arial','fontsize',fs)
    
    set(gca,'YTick',0.75:0.25:1.25)
    if (mod(ii,4)==1)==1
        set(gca,'YTickLabel',0.75:0.25:1.25, 'FontName','Arial','fontsize',fs)
    else
        
        set(gca,'YTickLabel',[]);
    end
    if ii==1
        text(0.1,0.83,'regular','units','normalized','color',col{1},'FontName','Arial','fontsize',fs)
        text(0.1,0.72,'permuted','units','normalized','color',col{2},'FontName','Arial','fontsize',fs)
    end
    if ii>2
        text(0.3,0.7,'bursty','color','red','units','normalized','FontName','Arial','fontsize',fs)
    end
    op=get(gca,'OuterPosition');
    set(gca,'OuterPosition',[op(1) op(2)+0.05 op(3) op(4)-0.05]); % OuterPosition = [left bottom width height]
    
    set(gca,'LineWidth',0.7,'TickLength',[0.025 0.025]);
end

axes
%text(-0.15,1.05,'A','units','normalized','FontName','Arial','fontsize',fs,'fontweight','bold')
h1 = xlabel ('frequency (Hz)','units','normalized','Position',[0.5,-0.04,0],'FontName','Arial','fontsize',fs);
h2 = ylabel ('power spectrum','units','normalized','Position',[-0.08,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h1,'visible','on')
set(h2,'visible','on')

%set(gcf, 'Position', get(0, 'Screensize'));
set(H, 'Units','centimeters', 'Position', pos_vec2) % size of the figure                                                     % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec2(3), pos_vec2(4)])      % for saving in the right size

if savefig2==1
    saveas(H,[savefile,figname2],'pdf');
end

%% plot ps and lower+upper bound of poisson for a specified session, 8 neurons
%{
cs=5;                                  % choose a session
N=size(ps_sess{cs},1);
vec=3:nf;
x=f(vec);
pos_vec2=[0,0,18,14];

H=figure('name',figname2,'visible','on');
for ii=1:8
    
    z1=lower{cs}(ii,vec);
    z2=upper{cs}(ii,vec);
   
    subplot(2,4,ii)
    hold on
    
    patch([x fliplr(x)], [z1 fliplr(z2)], col{2},'FaceAlpha',0.3,'EdgeColor',col{2},'EdgeAlpha',0.3,'Linewidth',0.5)
    plot(x,ps_sess{cs}(ii,vec),col{1},'linewidth',lw-0.5)
    hold off
    
    xlim([0,500])
    ylim([0.6,1.4])
    
    grid on
    
    set(gca,'XTick',0:250:500)
    if ii<N-4
        set(gca,'XTickLabel',[])
    else
        set(gca,'XTickLabel',0:250:500, 'FontName','Arial','fontsize',fs)
    end
    
    set(gca,'YTick',0.75:0.25:1.25)
    if (mod(ii,4)==1)==1
        set(gca,'YTickLabel',0.75:0.25:1.25, 'FontName','Arial','fontsize',fs)
    else
        
        set(gca,'YTickLabel',[]);
    end
    if ii==1
        text(0.1,0.83,'regular','units','normalized','color',col{1},'FontName','Arial','fontsize',fs)
        text(0.1,0.72,'permuted','units','normalized','color',col{2},'FontName','Arial','fontsize',fs)
    end
    
    set(gca,'LineWidth',1.0,'TickLength',[0.025 0.025]);
end

axes
%text(-0.15,1.05,'A','units','normalized','FontName','Arial','fontsize',fs,'fontweight','bold')
h1 = xlabel ('frequency (Hz)','units','normalized','Position',[0.5,-0.08,0],'FontName','Arial','fontsize',fs);
h2 = ylabel ('power spectrum','units','normalized','Position',[-0.08,0.5,0],'FontName','Arial','fontsize',fs);
set(gca,'Visible','off')
set(h1,'visible','on')
set(h2,'visible','on')

%set(gcf, 'Position', get(0, 'Screensize'));
set(H, 'Units','centimeters', 'Position', pos_vec2) % size of the figure                                                     % size of the figure
set(H,'PaperPositionMode','Auto','PaperUnits', 'centimeters','PaperSize',[pos_vec2(3), pos_vec2(4)])      % for saving in the right size

if savefig2==1
    saveas(H,[savefile,figname2],'pdf');
end
%}
