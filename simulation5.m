clear all
close all
clc
warning off  
echo off
format short g
% Assume minimum distance in X, and Y between
% stations

string_table=[];   

% Assume number of stations n

ppmatrix=[];   

n=20;
m=30;

% All the nodes with zero power
% pick any one to consider it max
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax=randi(n);   

load fixedTopology_n20_m30_lim15
neighbor_limit=25;  
dij=dist(cen1');   
% assuming a value to alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=2;    
         
Band=2000;             

lamdamV=.06*Band:.06*Band;   

% Threshold values
%%%%%%%%%%%%%%%%%%%
% Thr=(0.32*m)^2:(0.05*m)^2:(0.9*m)^2;  
Thr=0:100:300;     
Thr_count=0;

for vv=1:3   

Thr_count=Thr_count+1;

countLamda=0;  
counter=0;       

for lamdam=.02*Band:.02*Band          

lamdam    
countLamda=countLamda+1;  
    
% the power matrix
%%%%%%%%%%%%%%%%%%%%%%
pij=zeros(n);   

for n_runs=1:10                  

% cen1=stationsN(n,m,neighbor_limit);   


% Generate random numbers with mean lamdam and 
% variance lamdam/2
% Pick a request lamda_s_d
% 
%   
% Calculate the stations within range to each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k=1:n 
    ss=['station',num2str(k)];   
    TH=[];
    P1=cen1(k,:);
    neighbor=find( (dij(k,:) > 0) & (dij(k,:) <= neighbor_limit));  
    if (length(neighbor) >=1)
        for kk=1:length(neighbor)
            P2=cen1(neighbor(kk),:);   
            theta=theta1(P1,P2);
            TH(kk)=theta;  
        end;
        eval([ss,'.neighborsNo = neighbor;']);   
        eval([ss,'.neighborsAng = TH;']);  
    else
        eval([ss,'.neighborsNo = [];']);   
        eval([ss,'.neighborsAng = [];']);  
    end;    
end;  

% Calculate the stations within range to each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xij=zeros(n);
for k=1:n 
    ss=['station',num2str(k)];   
    k1=eval([ss,'.neighborsNo']);   
    xij(k,k1)=1;
    node_degree(k)=length(eval(['station',num2str(k), ...
        '.neighborsNo;']));        
end;       
clear ss dd rr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed Number pf requests for all runs
% for threshold analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Req_Bank_n_20_B_500_small_lamdam;   
eval(['reqmatAll=reqmat_run_',num2str(n_runs),';']);     
Ind1=find(reqmatAll(:,4)==lamdam);  
reqmat=reqmatAll(Ind1,:);    
reqmat=reqmat(:,1:3);   
clear reqmat_run_* reqmatAll   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed Number pf requests for all runs
% for threshold analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initilization
%%%%%%%%%%%%%%%%

consP=0;   
Tlamda_sd=0;       


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% consider number of requests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(['status_Thr_',num2str(Thr_count),'_run_',num2str(n_runs), ...
    '_lamda_',num2str(countLamda),'=[];']);       
eval(['string_table_Thr_',num2str(Thr_count),'_run_', ...
        num2str(n_runs),'_lamda_',num2str(countLamda),'=[];']);       

for nr=1:size(reqmat,1)             
    
% Pick any request    
%%%%%%%%%%%%%%%%%%%%
% kk=randi(size(reqmat,1));   

lamda_sd=reqmat(nr,1);   
s=reqmat(nr,2);   
d=reqmat(nr,3);   

% Generate the power vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PVij=ppij';    
PVij=transpose(dij.^alpha);    
Vxij=xij';   
PVij=PVij(:);   
Vxij=Vxij(:);   
PVij=PVij';    
Vxij=Vxij';     

CC=PVij;     
A_L_eq=[];   
B_L_eq=[];    
A_eq=[];
B_eq=[];   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The maximum power level Pmax
% Or use the total power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power constraint
% The total power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Power constraint as in the old paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TPower=sum(sum(triu(pij).*triu(xij)));    
% P>=pi>=dij^alpha*xij   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPmax=3000;   
ppi=sum(pij');    
if max(ppi) > 0   
    max_node=find(ppi==max(ppi));   
else
    max_node=randi(n,1);
end;       
max_node=max_node(1);   
% P_constraint=ppi(max_node)-Thr(n_runs) 
if max(ppi)>0
    TThr=Thr(vv)
    P_constraint=mean(ppi(find(ppi)))+Thr(vv)                    
else
    P_constraint=Thr(vv)                
end;   
for i=1:n    
    % pci is the consumed power in node i   
    pci=[zeros(length(1:(i-1)),n);pij(i,:)+dij(i,:).^alpha; ...
            zeros(length((i+1):n),n)];     
    % For the consumed power at node i
    pci=pci';
    pci=pci(:); 
    pci=pci';     
    % Variance constraint OUR invention
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % To reduce the variance
    % pi>PPmin           
    cond1=(ppi(s) < P_constraint ) && (max(ppi) > 0);   
    if  cond1 
        A_L_eq=[A_L_eq;pci];          
        B_L_eq=[B_L_eq;P_constraint];    
    end;       
    for j=1:n
        A_L_eq=[A_L_eq;pci];          
        ddij=PPmax-(dij(i,j).^alpha)*xij(i,j);   
        B_L_eq=[B_L_eq;ddij];    
        A_L_eq=[A_L_eq;pci];          
        B_L_eq=[B_L_eq;PPmax];                    
    end   
end;    
A_L_eq1=A_L_eq;    
B_L_eq1=B_L_eq; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delay constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta_sd   
% less than or equal constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume a value for delta_sd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_sd=ceil(2*n/3);  
A_L_eq=[A_L_eq;ones(1,n*n)];    
B_L_eq=[B_L_eq;delta_sd];     

A_L_eq2=ones(1,n*n);    
B_L_eq2=delta_sd; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandwidth constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Bandwidth B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AAA=[];   
for i=1:n   
    AA=zeros(n);  
    AA(i,:)=1;
    AA(:,i)=AA(:,i)+1;
    AA=AA';
    AA=AA(:);   
    AA=AA';
    AAA=[AAA;AA];   
end;    
A_L_eq=[A_L_eq;AAA];    
B_L_eq=[B_L_eq;(Band/lamda_sd)*ones(n,1)];    

A_L_eq3=AAA;    
B_L_eq3=(Band/lamda_sd)*ones(n,1);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New set of constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA1=[];   
for i=1:n   
    AA=zeros(n);  
    AA(i,:)=1;
    AA(:,i)=AA(:,i)-1;
    AA=AA';
    AA=AA(:);   
    AA=AA';
    AA1=[AA1;AA];   
    if s==i
        BB=1;   
    elseif d==i   
        BB=-1;  
    else
        BB=0;   
    end;      
    B_eq=[B_eq;BB];         
end;    
A_eq=[A_eq;AA1];    


% Another set of less than or equal   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AA=eye(n*n);    
B=xij';
B=B(:);
A_L_eq=[A_L_eq;AA];    
B_L_eq=[B_L_eq;B];     

A_L_eq4=AA;    
B_L_eq4=B;     

% Mixed integer linear programming function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb=zeros(1,n*n);   
ub=ones(1,n*n);   
e=2^-24;   
M=[1:n*n];   
t1=clock;   
t1=t1(4)*60*60+t1(5)*60+t1(6);   
save tt1 t1;   

clear  AA AA1 AAA B  BB  

[x,minpower,status]=IP1(CC,A_L_eq,B_L_eq,A_eq,B_eq,lb,ub,M,e); 

if status ==1   
    xres=reshape(x,n,n);   
    xres=xres';   
    last=s;   
    pp1=[last]; 
    while (last ~= d)
        last=find(xres(last,:) == 1);
        pp1=[pp1,last];
    end;    
    s
    d
    pp1    
% Reducing the limits for each tranmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kp=1:length(pp1)-1 
    ks=pp1(kp);  
    kd=pp1(kp+1); 
    ss=['station',num2str(ks)];   
    TH=[];
    P1=cen1(ks,:);
    shr=abs(sqrt(pij(ks,kd))*0.6);
    neighbor=find( (dij(ks,:) > 0) & (dij(ks,:) <= neighbor_limit-shr));  
    if (length(neighbor) >=1)
        for kk=1:length(neighbor)
            P2=cen1(neighbor(kk),:);   
            theta=theta1(P1,P2);
            TH(kk)=theta;  
        end;
        eval([ss,'.neighborsNo = neighbor;']);   
        eval([ss,'.neighborsAng = TH;']);  
    else
        eval([ss,'.neighborsNo = [];']);   
        eval([ss,'.neighborsAng = [];']);  
    end;    
end;  

% Calculate the stations within range to each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xij=zeros(n);
for k=1:n 
    ss=['station',num2str(k)];   
    k1=eval([ss,'.neighborsNo']);   
    xij(k,k1)=1;
    node_degree(k)=length(eval(['station',num2str(k), ...
        '.neighborsNo;']));        
end;       
    
end;    
St1=['Lamdam_',num2str(countLamda)];   
lamdam
status
node_Load=zeros(1,n);   
if status ==1  
    eval(['status_Thr_',num2str(Thr_count),'_run_',num2str(n_runs), ...
    '_lamda_',num2str(countLamda),'=[status_Thr_',num2str(Thr_count), ...
    '_run_',num2str(n_runs),'_lamda_',num2str(countLamda),',1];']);       
    path_string=[];  
    for k=1:length(pp1)-1
        path_string=[path_string,num2str(pp1(k)),'-->' ];
    end;
    path_string=[path_string,num2str(pp1(end))];    
    counter=counter+1;   
    consP=consP+minpower;  
    pij=pij+xres.*(dij.^alpha);  
    Tlamda_sd=Tlamda_sd+lamda_sd; 
    node_Load(pp1)=node_Load(pp1)+lamda_sd;   
    eval(['table1.lamdam(',num2str(counter),')=', ...   
        num2str(lamdam),';']);     
    eval(['table1.power(',num2str(counter),')=',num2str(minpower),';']);     
    eval('table1.Loads(counter,:)=node_Load;');     
    eval('table1.Degree(counter,:)=node_degree;');        
    eval(['string_table_Thr_',num2str(Thr_count),'_run_', ...
        num2str(n_runs),'_lamda_',num2str(countLamda), ...
        '=strvcat(string_table_Thr_',num2str(Thr_count),'_run_', ...
        num2str(n_runs),'_lamda_',num2str(countLamda),',path_string);']);           
%    figstr=['fig',num2str(counter),'_s_',num2str(s),'_d_',num2str(d)];  
%    saveas(gcf,figstr,'bmp')    
    ppi=sum(pij');   
%     imax=find(ppi == max(ppi));   
%     imax=imax(1);   
    eval(['pmatrix_run_',num2str(n_runs),'_lamda_', ...
        num2str(countLamda),'_req_',num2str(nr),'.power=ppi']);   
else 
    eval(['status_Thr_',num2str(Thr_count),'_run_',num2str(n_runs), ...
    '_lamda_',num2str(countLamda),'=[status_Thr_',num2str(Thr_count), ...
    '_run_',num2str(n_runs),'_lamda_',num2str(countLamda),',0];']);   
    eval(['string_table_Thr_',num2str(Thr_count),'_run_', ...
        num2str(n_runs),'_lamda_',num2str(countLamda), ...
        '=strvcat(string_table_Thr_',num2str(Thr_count),'_run_', ...
        num2str(n_runs),'_lamda_',num2str(countLamda),',[]);']);             
end;    


end   % lamda_s_d   (request)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
lamdaCell={'0.05B','0.1B'};

for k=1:length(lamdamV)
    ind=find(table1.lamdam==lamdamV(k));   
    if ~isempty(ind)    
        Loads=table1.Loads(ind,:);   
        Loads=sum(Loads);  
        max_LoadsALL(n_runs,k)=max(Loads);   
        min_LoadsALL(n_runs,k)=min(Loads);   
        avg_LoadsALL(n_runs,k)=mean(Loads);   
        Degrees=table1.Degree(ind,:);   
        Degrees=Degrees(end,:);  
        max_DegreesALL(n_runs,k)=max(Degrees);   
        min_DegreesALL(n_runs,k)=min(Degrees);   
        avg_DegreesALL(n_runs,k)=mean(Degrees);   
        Powers=table1.power(ind);  
        max_PowerALL(n_runs,k)=max(Powers);   
        min_PowerALL(n_runs,k)=min(Powers);   
        avg_PowerALL(n_runs,k)=mean(Powers);   
    end;
end; 
ppi
end % n_runs  
Thr(vv)
Var1(countLamda,vv)=var(ppi/max(ppi))         

end   % lamda_m

end;  % vv end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Arranging the requests
%%%%%%%%%%%%%%%%%%%%%%%%%

for k1=1:n_runs
    for k2=1:length(lamdamV)
        ss=['pmatrix_run_',num2str(k1),'_lamda_', ...
            num2str(k2),'_req_'];
        sss1=whos([ss,'*']);   
        sss2=[];
        for k3=1:size(sss1,1)
            sss2=strvcat(sss2,sss1(k3).name);
        end;   
        num_s=[];   
        for k3=1:size(sss1,1)
            num_s(k3)=str2num(sss2(k3,length(ss)+1:end));
        end        
        num_s=sort(num_s);   
        sst=[];  
        for k3=1:size(sss1,1)
            sst=strvcat(sst,['pmatrix_run_',num2str(k1),'_lamda_', ...
                num2str(k2),'_req_',num2str(num_s(k3))]);  
        end;
        eval(['run_',num2str(k1),'_lamda_',num2str(k2),'.Names=sst;']);                
    end;  
end;
  

for k3=1:countLamda
for k1=1:length(Thr)  
    tthr1=0;  
    tthr2=0;  
    tthr3=0;      
    for k2=1:n_runs           
        vv1=eval(['status_Thr_',num2str(k1), ...
            '_run_',num2str(k2),'_lamda_',num2str(k3),';']);   
        tthr1=tthr1+length(find(vv1==0));   
        tthr2=tthr2+length(find(vv1==1));           
        tthr3=tthr3+length(vv1);           
    end;  
    LostV(k1)=tthr1;   
    SendV(k1)=tthr2;   
    PacketV(k1)=tthr3;       
end;        
LLostV(k3,:)=LostV;   
SSendV(k3,:)=SendV;   
PPacketV(k3,:)=PacketV;       
end;  



% Table Generation
%%%%%%%%%%%%%%%%%%%%%
for k3=1:countLamda   
    for k1=1:length(Thr)
        for k2=1:n_runs   
            eval(['sstatus=status_Thr_',num2str(k1),'_run_', ...
                num2str(k2),'_lamda_',num2str(k3),';']);    
            load Req_Bank_n_20_B_500_small_lamdam;   
            eval(['reqmatAll=reqmat_run_',num2str(k2),';']);     
            lamdam=lamdamV(k3);  
            Ind1=find(reqmatAll(:,4)==lamdam);  
            reqmat=reqmatAll(Ind1,:);    
            reqmat=reqmat(:,1:4);  
            clear reqmat_run_* reqmatAll;   
            rr1=num2str(reqmat);   
            Table1=[];  
            counter1=0;   
            for k4=1:length(sstatus)   
                if sstatus(k4) == 1   
                    counter1=counter1+1;   
                    eval(['Table1=strvcat(Table1,string_table_Thr_', ...
                        num2str(k1),'_run_',num2str(k2),'_lamda_', ...
                        num2str(k3),'(counter1,:));']);  
                else
                    eval(['Table1=strvcat(Table1,''  Lost '');']);
                end;
            end;  
            for k=1:size(Table1,1)
                TTable1(k,:)=[blanks(5),Table1(k,:)];    
            end;
            Table1=TTable1;   
            clear TTable1;   
            eval(['Table_Thr_',num2str(k1),'_run_',num2str(k2), ...
                '_lamda_',num2str(k3),'=[rr1,Table1];']);   
        end;
    end;
end;

% save sim_18_5_2010_ThresholdAnalysis_nr_constraints    




