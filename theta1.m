function Theta = theta1(P1, P2)
% Calculate the angle between two points P1 and P2

% Calculate the differences in coordinates
deltaX = P2(1) - P1(1);
deltaY = P2(2) - P1(2);

% Calculate the angle
if deltaX == 0
    if deltaY >= 0
        Theta = pi/2;
    else
        Theta = 3*pi/2;
    end
elseif deltaY == 0
    if deltaX >= 0
        Theta = 0;
    else
        Theta = pi;
    end
elseif deltaX > 0 && deltaY > 0
    Theta = atan(deltaY / deltaX);
elseif deltaX < 0 && deltaY > 0
    Theta = pi - atan(abs(deltaY) / abs(deltaX));
elseif deltaX < 0 && deltaY < 0
    Theta = pi + atan(abs(deltaY) / abs(deltaX));
else
    Theta = 2*pi - atan(abs(deltaY) / abs(deltaX));
end

% Convert angle to degrees
Theta = rad2deg(Theta);

end
