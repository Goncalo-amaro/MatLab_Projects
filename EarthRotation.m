function [X, Y, Z] = EarthRotation(X, Y, Z, t)

% Z(1,:) = -Z(2,:);
psi = -360*pi/180/24/60/60*t;
theta = 0;
phi = 0;
A = angle2dcm(psi, theta, phi,'ZXZ');
for j =1:1:length(X(1,:))
    for i = 1:1:length(X(1,:))
       R = A*[X(j,i);Y(j,i);Z(j,i)];
       X(j,i)=R(1);
       Y(j,i)=R(2);
       Z(j,i)=R(3);
    end
end

end

