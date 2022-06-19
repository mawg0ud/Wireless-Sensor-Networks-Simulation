clear all
close all
clc
format short g   
Band=500;   
counter=0;   
load centers_n_15_2;

string_table=[];   

neighbor_limit=80;  
Small_Stations(cen1,num1,limits,neighbor_limit);   
dij=zeros(n,n);   

for k=1:n   
    for kk=1:n   
        dij(k,kk)=sqrt((cen1(k,1)-cen1(kk,1))^2+(cen1(k,2)-cen1(kk,2))^2);    
    end;
end;

% Calculate the stations within range to each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
end;       



% assuming a value to alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=2;    
         
% Generate random numbers with mean lamdam and 
% variance lamdam/2
% Pick a request lamda_s_d
% 
%

lamdam=150;  
rr = lamdam + sqrt(lamdam/2).*randn(10000,1);

% var(r),mean(r)

% The traffic demand between node pair (s,d) lamda_sd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random request generation (s is row, d is column)
sd=randi([0 1],n);   
for k=1:n   
    sd(k,k)=0;
end;      

% Find the request in vector format ss and dd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dd,ss]=find(sd');   

while length(dd)<15
sd=randi([0 1],n);   
for k=1:n   
    sd(k,k)=0;
end;      

% Find the request in vector format ss and dd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dd,ss]=find(sd');   
end;    

% Select lamda_s_d
% Generate request matrix  
index1=1:length(dd);  
for k=1:length(dd)-3
    index=randi(length(rr));
    reqmat(k,1)=rr(index);   
    index=randi(length(dd));   
    reqmat(k,2)=ss(index);       
    reqmat(k,3)=dd(index);   
    index1=find(index1 ~= index);
    dd=dd(index1);       
    ss=ss(index1);    
    index1=1:length(dd);      
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Band=500:500:1500     

% the power matrix
%%%%%%%%%%%%%%%%%%%%%%
pij=dij.^alpha;    
    
% Initilization
%%%%%%%%%%%%%%%%

consP=0;   
Tlamda_sd=0;       

%%%%%%%%%%%%%%%%        
    
% consider number of requests=6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nr=1:12     
    
% Pick any request    
%%%%%%%%%%%%%%%%%%%%
% kk=randi(size(reqmat,1));   

lamda_sd=reqmat(nr,1);   
s=reqmat(nr,2);   
d=reqmat(nr,3);   

% Generate the power vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%
PVij=pij';   
PVij=PVij(:); 
PVij=PVij';    

% Assume that we have one request s -->d   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The coefficient of the optimized function
% is PVij which is C .....
% We use the total power instead
% The unknowns are 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

TPower=sum(sum(triu(pij).*triu(xij)));    
A_L_eq=[A_L_eq;CC];    
B_L_eq=[B_L_eq;TPower];     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delay constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta_sd   
% less than or equal constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume a value for delta_sd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_sd=4;  
A_L_eq=[A_L_eq;ones(1,n*n)];    
B_L_eq=[B_L_eq;delta_sd];     

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
B_L_eq=[B_L_eq;Band/lamda_sd*ones(n,1)];    


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

% Mixed integer linear programming function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb=zeros(1,n*n);   
ub=ones(1,n*n);   
e=2^-24;   
M=[1:n*n];   
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
end;    

if status ==1  
    path_string=[];  
    for k=1:length(pp1)-1
        path_string=[path_string,num2str(pp1(k)),'-->' ];
    end;
    path_string=[path_string,num2str(pp1(end))];    
    counter=counter+1;   
    consP=consP+minpower;  
    pij=pij+xres.*pij;  
    Tlamda_sd=Tlamda_sd+lamda_sd;       
    table1(counter,1)=n;   
    table1(counter,2)=s;   
    table1(counter,3)=d;   
    table1(counter,4)=Tlamda_sd;   
    table1(counter,5)=Band;   
    table1(counter,6)=consP;       
    string_table=strvcat(string_table,path_string);           
    figstr=['fig',num2str(counter),'_s_',num2str(s),'_d_',num2str(d)];  
    saveas(gcf,figstr,'jpg')    
end;    

end   % lamda_s_d   

end   % Bandwidth    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the data to a text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('res_simulation3.txt','w');   
fprintf(fid,'%70s\n\n\n',datestr(now,0));   

fprintf(fid,'\t\t Ad Hoc, Power Saving Analysis \n');  

fprintf(fid,'=================== Data Analysis ====================\n');  
fprintf(fid,'==========================================================\n\n'); 
fprintf(fid,'n\t s\t d\t cum. lamda\t Bandwidth \t\t  route\t\t\t\t\t Power\r\n');   
fprintf(fid,'--- --- --- ----------\t -------- -----------------------\t-------------\t \r\n');      
for k=1:size(table1,1)  
   fprintf(fid,'%2d  %2d  %2d  %10.4f  %7d \t  %20s  %10d \r\n', ...
       table1(k,1),table1(k,2),table1(k,3),table1(k,4),table1(k,5), ...
            string_table(k,:),table1(k,6));  
end;     

fclose(fid);     


save route1    

