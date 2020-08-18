% compute distance of weights w S+C and w C 

clear all 
close all
clc 

saveres=0;
type=3;

Krange=[300,400,500];

%%
namet={'sign','info','bursty'};
namei={'S+C','C'};
variable={'plus_alltr','tag_info','bursty'};

task=['distance weights S+C/C for ',namet{type}];
disp(task)

%% load results

addpath('/home/veronika/synced/transfer_result/weight/weight_comp2/')
addpath('/home/veronika/synced/transfer_result/tags/')

loadname=['tag_',namet{type},'_V1_400'];
load(loadname)

%%

D=cell(3,2);

for int=1:3
    
    K=Krange(int);
    
    clear w_comp
    loadname=['weight_comp_',sprintf('%1.0i',K)];
    load(loadname)
    w_abs=cellfun(@(x) abs(x), w_comp, 'UniformOutput', false);
    
    nbses=size(w_abs,2);
    ncv=size(w_abs{1},1);
    %% collect weights across sessions
    
    w1=[];
    w2=[];
    for sess=1:nbses
        x=w_abs{1,sess};
        w1=cat(2,w1,x);
        
        y=w_abs{2,sess};
        w2=cat(2,w2,y);
    end
    
    %% compute distance between weight d(w_n^{s+c},w_n^{c}) for each of the two groups
   
    tags=cell(2,1);
    tags{1}=cell2mat(eval(variable{type}));               % bursty
    tags{2}=cell2mat(eval(variable{type}))==0;             % non-bursty
    %%
    for u=1:2
        
        idx=find(tags{u}'==1);              % index of neurons from a group
        Nuse=length(idx);                   % number of neurons in a  group
        
        d=zeros(ncv,Nuse); 
        w1u=w1(:,idx);
        w2u=w2(:,idx);
        
        for n=1:Nuse
            d(:,n)=abs(w1u(:,n)-w2u(:,n));
        end
        
        D{int,u}=mean(d); % average across cv
    end
end
%%
if saveres==1
    
    savename=['w_distance_',namet{type}];
    savefile='/home/veronika/synced/transfer_result/weight/w_distance/';
    save([savefile,savename],'D','Krange')
    %clear all
end    
