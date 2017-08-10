%% Homework 4 Exercise 1.a.)

xvec = 0:.05:.4; yvec = 0:.05:.4; zvec = 0:.05:.4;

%called it
c = [.1; 1.8; 11.8];
[x0, y0, z0] = meshgrid(xvec + c(1), yvec + c(2), zvec+ c(3));
%this sort of spits out coordinate reprsentations of planes, so that we
%then have coordinates at every xyz point in a cube in increments of 2.

% Now creating an array of initial conditions where the first array is my x
% coordinates etc. So this is actually an indexed array of all the meshgrid
% coords. I guess geometrically it's a tenser, but computationally it's
% just a list of 3d matrix values.
% yIC(1, :, :, :) = x0;
% yIC(2, :, :, :) = y0;
% yIC(3, :, :, :) = z0;

x_v = x0(:); y_v = y0(:);z_v = z0(:);
A1 = [transpose(x_v);transpose(y_v);transpose(z_v)];

% plot3(x_v, y_v, z_v, '.')

%% Exercise 1.b.)
dt = .01;
t = [0:dt:7];
list_coords = zeros(701, 3, 729);
list_coords(1,:,:) = A1;
y_in = squeeze(list_coords(1,:,:));
counter = 0;
time = 0;
for i = 1:t(end)/dt %this gives the number of steps in tspan
    time = i * dt;
    y_out = rk4_singleStep(@(t,y) chaos_ode(t, y), dt, time, y_in);
    list_coords(i+1,:,:) = y_out;
    y_in = y_out;
%     plot3(y_out(1,:), y_out(2,:), y_out(3,:), 'r.' , 'LineWidth', 2, 'Markersize', 10)
%     drawnow
    counter = counter + 1;
end

A2 = squeeze(list_coords(:,1,:));
A3 = squeeze(list_coords(:,2,:));
A4 = squeeze(list_coords(:,3,:));

% figure
% plot3(A2(end,:),A3(end,:),A4(end,:));
% plot3(A2(:,1:20),A3(:,1:20),A4(:,1:20));
% for i=1:length(t)
%     plot3(A2(i,:),A3(i,:),A4(i,:),'.');
%     view(20+360*i/(7/.01), 40);
%     axis([-.5,.5,-10,10,0,50])
%     drawnow
%     pause(dt)
% end