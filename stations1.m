% Demonstrate the random spread of a group of stations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
format short g
% Assume minimum distance in X, and Y between
% stations

string_table=[];   
% Assume number of stations n
counter=0; 
tt=1;  
while tt ==1   
n=input(['Enter the number of nodes you assumed', ...
    'the value should be integer more than 1 \n\n']);    
n=round(n);
if n>1
tt=0;
end;
end;

        
xx=10:500:18000;   
yy=10:500:18000;   
xmin=xx(1)-500;  
xmax=xx(end)+500;  
ymin=yy(1)-500;  
ymax=yy(end)+500;  
limits=[xmin xmax ymin ymax];  
% Random locations selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:n
    xr=randi(length(xx));   
    xc(k)=xx(xr);  
    xi=find(xx ~= xc(k));
    xx=xx(xi);   
    yr=randi(length(xx));   
    yc(k)=yy(yr);  
    yi=find(yy ~= yc(k));
    yy=yy(yi);
end;    
    
cen1=[xc',yc'];   
num1=1:n;
neighbor_limit=14000;  
stations(cen1,num1,limits,neighbor_limit);  

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

