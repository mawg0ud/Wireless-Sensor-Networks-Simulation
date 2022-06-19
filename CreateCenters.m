% Demonstrate the random spread of a group of stations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
format short g
% Assume minimum distance in X, and Y between
% stations
t1=clock;   
string_table=[];   
% Assume number of stations n
counter=0;   
step1=40;  
n=25;    
xx=10:step1:18000;   
yy=10:step1:18000;   
xmin=xx(1)-20*step1;  
xmax=xx(end)+20*step1;  
ymin=yy(1)-20*step1;  
ymax=yy(end)+20*step1;  
limits=[xmin xmax ymin ymax];  
% Random locations selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:n
    xr=randi(length(xx));   
    xc(k)=xx(xr);     
    xstep=xc(k)-10*step1:step1:xc(k)+10*step1;   
    for kk=1:length(xstep)
        xi=find(xx ~= xstep(kk)); 
        xsize=size(xx)
        xx=xx(xi);   
    end;      
    yr=randi(length(yy));   
    yc(k)=yy(yr);     
    ystep=yc(k)-10*step1:step1:yc(k)+10*step1;   
    for kk=1:length(ystep)
        yi=find(yy ~= ystep(kk));   
        ysize=size(yy)        
        yy=yy(yi);   
    end;      
end;    
    
cen1=[xc',yc'];   
num1=1:n;
neighbor_limit=18000;  
stations(cen1,num1,limits,neighbor_limit);  
grid
t2=clock; 
t3=t2-t1  




% save centers_n_25_6