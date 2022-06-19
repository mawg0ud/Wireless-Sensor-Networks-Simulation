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
