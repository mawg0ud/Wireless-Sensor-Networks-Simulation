function [Theta]=theta1(P1,P2)
% function [theta]=theta1(P1,P2)
% This function calculate the angle 
% between two station P1 sender, P2
% Receiver
deltaX=P2(1)-P1(1);  
deltaY=P2(2)-P1(2);  

if (deltaX == 0) && (deltaY >= 0)
    Theta = pi/2; 
elseif (deltaX == 0) && (deltaY < 0)
    Theta=3*pi/2;  
elseif (deltaX >= 0) && (deltaY == 0)
    Theta=0;      
elseif (deltaX < 0) && (deltaY == 0)
    Theta=pi;      
elseif (deltaX > 0) && (deltaY > 0)
    Theta=atan(abs(deltaY/deltaX)); 
elseif (deltaX < 0) && (deltaY > 0) 
    Theta=pi-atan(abs(deltaY/deltaX));     
elseif (deltaX < 0) && (deltaY < 0) 
    Theta=pi+atan(abs(deltaY/deltaX));         
else (deltaX > 0) && (deltaY < 0); 
    Theta=2*pi-atan(abs(deltaY/deltaX));             
end;      
Theta=Theta*180/pi;  