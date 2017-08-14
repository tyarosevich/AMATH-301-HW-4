%% Exercise 4b. a.)
clc;clear all;close all;
load imag_data.mat

A1 = sum(B, 2)/40;

%imshow(uint8(reshape(A1,200,175)))

%% Exercise 4b. b.)

A = B;
for j = 1:40
    A(:,j) = A(:,j) - A1;
end

[U, S, V] = svd(A, 'econ');
s = svd(A,'econ');
A2 = s(1:10);

%% Exercise 4b. c.)

%plot(1:40, s)

x = transpose(U(:,1)) * A;
y = transpose(U(:,2)) * A;
z = transpose(U(:,3)) * A;
A3 = x;

%% Exercise 4b. d.)

% plot3(x(:,1:20), y(:,1:20),z(:,1:20), 'bo')
% hold on
% plot3(x(:,21:40), y(:,21:40),z(:,21:40), 'ro')

%% Exercise 4b. e.)

u_2 = u - A1;

A_final = [A u_2];
x = transpose(U(:,1)) * A_final;
y = transpose(U(:,2)) * A_final;
z = transpose(U(:,3)) * A_final;

% plot3(x(:,1:20), y(:,1:20),z(:,1:20), 'bo')
% hold on
% plot3(x(:,21:40), y(:,21:40),z(:,21:40), 'ro')
% plot3(x(:,41), y(:,41),z(:,41), 'gs')

%imshow(uint8(reshape(u,200,175)))





