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
%% Exercise 2.a.)
clc;clear all;close all;

f = @(x,r) r.*x.*(1-x);
A5 = zeros(31, 4);
rList = [2.5; 3.2; 3.52; 4];
A5(1,:) = .7;
for j  = 1:4
    r = rList(j);
    for i = 2:31
        A5(i, j) = f(A5(i-1, j),r);
    end
end
% t = 0:30;
% plot(t, A5(:,3))
% axis([-5 40 0.25 1])
% grid on
% grid minor

P_r1 = 1; P_r2 = 2;P_r3 = 4;
A6 = [P_r1;P_r2;P_r3];
%% Exercise 2.b.)
r = 2.5:.001:4;
A7 = zeros(501, 1501);
A7(1,:) = .7;
for i = 2:501
    A7(i,:) = f(A7(i-1,:),r);

end

A7_final = A7(402:end,:);
%plot(r,A7_final,'k.','Markersize',1);

%% Exercise 3.a.)

dt = .004;dx = .01;phi = .1;L = 1;
steps = L/dx - 1;
A = (1 - 2*((phi^2 * dt)/(dx^2)));
B = (phi^2 * dt)/(dx^2);
M = diag(zeros(1,steps)+ A) + diag((zeros(1,steps - 1)+B), -1) + diag((zeros(1,steps -1) +B), 1);



%% Exercise 3.b.)

M_eigs = eig(M);
A11 = max(abs(M_eigs));

%% Exercise 3.c.)

x_i = dx:dx:L-dx;
u_0 = exp(-200*(x_i -.5).^2);

A12 = transpose(u_0);

%% Exercise 3.d.)

U_List = zeros(99,501);
U_List(:,1) = A12;
for n = 1:2/dt
    t = n*dt;
    U_List(:, n+1) = M*U_List(:,n);
end

u_bounds = zeros(1, 501);

U_List = [u_bounds;U_List;u_bounds];

A13 = U_List;

%% Exercise 3.e.)

% for j=1:(2/dt+1)
% plot(linspace(0,L,L/dx+1)',A13(:,j));
% axis([0,L,0,1]);
% pause(dt);
% end
%     
    










