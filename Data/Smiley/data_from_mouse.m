% data_from_mouse.m
% MATLAB code to take graphical input from the mouse, plot, and write to file


figure
n = 100;
coordinates = zeros(n,2);
hold on
for i=1:n
[x, y] = ginput(1);
coordinates(i,:) = [x, y];
plot(coordinates(:,1), coordinates(:,2), '+');
end
hold off
csvwrite('smiley.dat',coordinates)
