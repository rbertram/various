% Fig6.m

% Used to make Figure 6A-D from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% This code solves a simple DDE system modelling mutual positive feedback of 2
% molecules. After, we plot phase plane trajectories for specific tau and initial conditions.
% X -> Y activation of Y by X
% Y -> X activation of X by Y

% Parameters_PosPos and SteadyState_Solutions_PosPos.m must be in current directory!!

% Delay Equations:
% dx/dt = (alphaX*ylag^n) / (1 + ylag^n) - x;
% dy/dt = (alphaY*xlag^n) / (1 + xlag^n) - y;

% Author: Christopher Ryzowicz
%-------------------------------------------------------------------------
clear all
clc
% close all

%%% Timespan and delay
tau = 6;% Delay value
Tstart = 0; % Start time of solution
Tend = 2000; % End time of solution

%%% Parameters
Parameters_PosPos

%%% Steady State Solutions
[ss_x, ss_y] = SteadyState_Solutions_PosPos();
ss_x = real(ss_x(imag(ss_x) == 0)); % Extract only the real steady states
ss_y = real(ss_y(imag(ss_y) == 0)); % Extract only the real steady states

%%% History constants (x_0,y_0)
hist = {[0.85 0.61], [0.865 0.625]};

%%% Generating points
x1 = 0:0.01:2.15;

%%% Plot the stable manifold
% !!!Stable manifold will have to be approximated by a straight line using the eigenvectors of Jacobian near the saddle point!!!
figure(1),
plot(x1,1.36-x1,'k--','LineWidth',2),hold on
figure(2),
plot(x1,1.36-x1,'k--','LineWidth',2),hold on

%%% Plot steady state points with filled circles signifying stable s.s and triangle as saddle point
figure(1),
plot(ss_x(1:2:3),ss_y(1:2:3),'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(ss_x(2),ss_y(2),'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
figure(2),
plot(ss_x(1:2:3),ss_y(1:2:3),'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(ss_x(2),ss_y(2),'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on


for i=1:length(hist)
    C = {['r','g']}; % Assign colors for trajectory plot

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol = dde23(@ddefunc, tau, hist{i}, [Tstart Tend], opts, q);

    %%% Plot trajectory
    figure(i),
    plot(sol.y(1,:),sol.y(2,:),'Color',C{1}(i)), hold on
    xlabel('x')
    ylabel('y')
    box off
    set(gca,'FontSize',24)
end

%%% Adjust axes
figure(1), axis([0 2 0 2])
xticks([0 1 2])
yticks([0 1 2])
plot(hist{1}(1),hist{1}(2),'pentagram','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
title('Figure 6C')
figure(2), axis([0 2 0 2])
xticks([0 1 2])
yticks([0 1 2])
plot(hist{2}(1),hist{2}(2),'pentagram','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
title('Figure 6D')

%---------------------------------------------------------------------------------------------
%%% Plot the stable manifold
% !!!Stable manifold will have to be approximated by a straight line using the eigenvectors of Jacobian near the saddle point!!!
figure(3),
plot(x1,1.36-x1,'k--','LineWidth',2),hold on
figure(4),
plot(x1,1.36-x1,'k--','LineWidth',2),hold on

%%% Plot steady state points with filled circles signifying stable s.s and triangle as saddle point
figure(3),
plot(ss_x(1:2:3),ss_y(1:2:3),'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(ss_x(2),ss_y(2),'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
figure(4),
plot(ss_x(1:2:3),ss_y(1:2:3),'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(ss_x(2),ss_y(2),'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on

delays = [4 6];
for j=1:length(delays)
    C = {['r','g']}; % Assign colors for trajectory plot

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol = dde23(@ddefunc, delays(j), [1.2 0.5], [Tstart Tend], opts, q);

    %%% Plot trajectory
    figure(j+2),
    plot(sol.y(1,:),sol.y(2,:),'Color',C{1}(j)), hold on
    xlabel('x')
    ylabel('y')
    box off
    set(gca,'FontSize',24)
end

%%% Adjust axes
figure(3), axis([0 2 0 2])
xticks([0 1 2])
yticks([0 1 2])
plot(1.2,0.5,'pentagram','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
title('Figure 6A')
figure(4), axis([0 2 0 2])
xticks([0 1 2])
yticks([0 1 2])
plot(1.2,0.5,'pentagram','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
title('Figure 6B')

%----------------------------------------------------------------------

function dout = ddefunc(t, out, yl, q)
%%% Define the variables
x = out(1);
y = out(2);
xlag = yl(1);
ylag = yl(2);

%%% Differential Equation
dx = (q.alphaX*ylag^q.n) / (1 + ylag^q.n) - x;
dy = (q.alphaY*xlag^q.n) / (1 + xlag^q.n) - y;

dout = [dx dy]';
end
