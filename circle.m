function [xc,yc]=circle(cen,radius)   
% function [xc,yc]=circle(cen,radius);   
% xc a vector of x coordinates ...
% yc a vector of y coordinates ...
% cen is the center of the circle
% radius of the circle
th1=0:.01:2*pi;   
xc=radius*cos(th1);   
yc=radius*sin(th1);   
xc = xc + cen(1);  
yc = yc + cen(2);  


