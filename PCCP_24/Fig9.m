% Fig9.m

% Used to make Figure 9B-C from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% DDE Formation from "Detection of multistability, bifurcations, and
% hysteresis in a large class of biological positive-feedback systems"
% Paper of the Cdc2-Wee1 protein interaction. Plotting phase-space and time series for ODE and DDE trajectory.

% Author: Christopher Ryzowicz
%-------------------------------------------------------------------------
clear all
clc
clf
%% FigureB-----------------------------------------------------------------
%%%Timespan and delay
tau = 8; % Delay value
Tstart = 0; % Start time of solution
Tend = 1000; % End time of solution

%%% History constants (x_0,y_0)
hist = {[0.25,0.62], [0.4,0.3]};

%%% Plot stable manifold (estimated via a straight line from the eigenvectors of Jacobian near saddle point)
x1 = 0:0.01:1.1;
figure(1),
plot(x1,0.1706+0.8877*x1,'k--','LineWidth',2),hold on

%%% Plot steady state points with filled circles signifying stable s.s and triangle as saddle point
figure(1),
plot(0.13573,0.99662,'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(0.9947,0.16816,'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(0.50577,0.61951,'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on

for i=1:length(hist)
    C = {['b','r']}; % Assign colors for trajectory plot

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol = dde23(@ddefunc, tau, hist{i}, [Tstart Tend], opts);

    %%% Solving ODE
    solODE = ode23s(@odefunc, [Tstart Tend], hist{i}, []);

    %%% Plot trajectory
    figure(1),
    plot(sol.y(1,:),sol.y(2,:),'LineWidth',1,'Color',C{1}(i)), hold on
    plot(solODE.y(1,:),solODE.y(2,:),'k','LineWidth',2), hold on
    xlabel('x')
    ylabel('y')
    box off
    set(gca,'FontSize',20)
end
figure(1), axis([0 1.1 0 1.1])
xticks([0 0.5 1])
yticks([0 0.5 1])

%% Figure 7C--------------------------------------------------------------

%%% Solving DDE
opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
sol2 = dde23(@ddefunc, tau, hist{2}, [Tstart Tend],opts);

%%% Plotting time course for Active Cdc2 and active Wee1 protein normalized time for number of oscillations by dividing by 2*tau
figure(2),
plot(sol2.x/(2*tau),sol2.y(1,:),'Color',[0 0.45 0.73],'LineWidth',2)
xlabel('oscillation number')
ylabel('x')
box off
set(gca,'FontSize',20)
axis([0 41 0 1])
xticks([0 10 20 30 40])
yticks([0 0.5 1])

%------------------------------------------------------------------------

function dout = ddefunc(t, y, yl)
%%% Differential equations function

% Define the variables
x1 = y(1);
y1 = y(2);
x1lag = yl(1,:);
y1lag = yl(2,:);

% Parameters from the Supplementary material of the paper
nu = 1;
alpha1 = 1;
alpha2 = 1;
beta1 = 200;
beta2 = 10;
gamma1 = 4;
gamma2 = 4;
k1 = 30;
k2 = 1;

% Differential Equation
dx1 = alpha1*(1 - x1) - (beta1*x1*(nu*y1lag)^gamma1) / (k1 + (nu*y1lag)^gamma1);
dy1 = alpha2*(1 - y1) - (beta2*y1*(x1lag)^gamma2) / (k2 + (x1lag)^gamma2);

dout = [dx1 dy1]';
end

%------------------------------------------------------------------------

function dout = odefunc(t, y)
%%% Differential equations function

% Define the variables
x1 = y(1);
y1 = y(2);

% Parameters from the Supplementary material of the paper
nu = 1;
alpha1 = 1;
alpha2 = 1;
beta1 = 200;
beta2 = 10;
gamma1 = 4;
gamma2 = 4;
k1 = 30;
k2 = 1;

% Differential Equation
dx1 = alpha1*(1 - x1) - (beta1*x1*(nu*y1)^gamma1) / (k1 + (nu*y1)^gamma1);
dy1 = alpha2*(1 - y1) - (beta2*y1*(x1)^gamma2) / (k2 + (x1)^gamma2);

dout = [dx1 dy1]';
end
