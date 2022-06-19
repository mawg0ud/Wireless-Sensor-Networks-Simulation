function []=Small_Stations(cen1,num1,limits,neighbor_limit)   
% function []=stations(cen1,num1,neighbor_limit)   
% cen1 is the center of the stations
% num1 is the numbers
figure('position',[100 100 850 600]);   
xxc=[];  
yyc=[];  
neighbor_xc=[];  
neighbor_yc=[];  
for k=1:length(num1)
    [xc yc]=circle([cen1(k,1) cen1(k,2)],7);    
    [xn yn]=circle([cen1(k,1) cen1(k,2)],neighbor_limit);        
    xxc=[xxc,xc'];   
    yyc=[yyc,yc'];   
    neighbor_xc=[neighbor_xc,xn'];  
    neighbor_yc=[neighbor_yc,yn'];      
end;    

xxc=[xxc,neighbor_xc];   
yyc=[yyc,neighbor_yc];   
plot(xxc,yyc),axis image,
axis([limits(1) limits(2) limits(3) limits(4)]);       
hold on
for k=1:length(num1)
    gt=text(cen1(k,1)-3,cen1(k,2)+1,num2str(num1(k)));   
    set(gt,'FontSize',12,'FontWeight','bold')
end;
xlabel('Units in meters');  
ylabel('Units in meters');  
