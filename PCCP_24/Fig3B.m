% Fig3B.m

% Used to make Figure 3E-F from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.


% This code solves a simple DDE system modelling mutual positive feedback of 2
% molecules. After, we plot 2 trajectories, one from ODE and other from
% DDE. In addition, we provide a sample time series in x for the DDE system.

% Parameters_PosPos.m and SteadyState_Solutions_PosPos.m must be in current directory!!

% X -> Y activation of Y by X
% Y -> X activation of X by Y

% Delay Equations:
% dx/dt = (alphaX*ylag^n) / (1 + ylag^n) - x;
% dy/dt = (alphaY*xlag^n) / (1 + xlag^n) - y;

% Author: Christopher Ryzowicz 
%-------------------------------------------------------------------------
clear all
clc
close all

%% Figure 3E----------------------------------------------------------------
%%% Timespan and delay
tau = 8;% Delay value
Tstart = 0; % Start time of solution
Tend = 1000; % End time of solution

%%% Parameters
Parameters_PosPos

%%% Steady State Solutions
[ss_x, ss_y] = SteadyState_Solutions_PosPos();
ss_x = real(ss_x(imag(ss_x) == 0)); % Extract only the real steady states
ss_y = real(ss_y(imag(ss_y) == 0)); % Extract only the real steady states

%%% History constants (x_0,y_0)
hist = {[0.2,1.1], [1.4,0.6]};

%%% Plot the stable manifold
%%% !!!Stable manifold will have to be approximated by a straight line using the eigenvectors of Jacobian near the saddle point!!!
x1 = 0:0.01:2.15;
figure(1),
plot(x1,1.36-x1,'k--','LineWidth',2),hold on

%%% Plot steady state points with filled circles signifying stable s.s and triangle as saddle point
figure(1),
plot(ss_x(1:2:3),ss_y(1:2:3),'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(ss_x(2),ss_y(2),'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on

for i=1:length(hist)
    C = {['r','b']}; % Assign colors for trajectory plot

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol = dde23(@ddefunc, tau, hist{i}, [Tstart Tend], opts, q);

    %%% Solving ODE
    solODE = ode23s(@odefunc, [Tstart Tend], hist{i}, [], q);

    %%% Plot trajectories
    figure(1),
    plot(sol.y(1,:),sol.y(2,:),'LineWidth',1,'Color',C{1}(i)), hold on
    plot(solODE.y(1,:),solODE.y(2,:),'k','LineWidth',2), hold on
    xlabel('x')
    ylabel('y')
    box off
    set(gca,'FontSize',20)
end

%%% Adjust axes
figure(1), axis([0 2 0 2])
xticks([0 1 2])
yticks([0 1 2])

%% Figure 3F-----------------------------------------------------------

%%% Solving DDE
opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
sol2 = dde23(@ddefunc, tau, hist{1}, [Tstart Tend],opts,q);

%------------------------------------------------------------------------

%%% Plot time series for DDE system
figure(2),
plot(sol2.x/(2*tau),sol2.y(1,:),'LineWidth',2)
xlabel('oscillation number')
ylabel('x')
box off
set(gca,'FontSize',20)
axis([0 45 0 1.6])
yticks([0 0.8 1.6])

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

function dout1 = odefunc(t, out, q)
%%% Define the variables
x = out(1);
y = out(2);

%%% Differential Equation
dx = (q.alphaX*y^q.n) / (1 + y^q.n) - x;
dy = (q.alphaY*x^q.n) / (1 + x^q.n) - y;

dout1 = [dx dy]';
end
