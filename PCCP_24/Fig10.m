% Fig10.m

% Used to make Figure 10B-C from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% DDE simulations of the gene/protein interaction toy model. We plot phase
% plane trajectories and time series for ODE and DDE models.
%
% Author: Christopher Ryzowicz
%-------------------------------------------------------------------------
clear all
clc
clf
%% Figure 8B--------------------------------------------------------------
%%%Timespan and delay
tau = 10; % Delay value
Tstart = 0; % Start time of solution
Tend = 1200; % End time of solution

%%% History constants (x_0,y_0)
% hist = {[0.2,0.6], [0.17,2.2]};
hist = {[0.22,1], [0.17,2.2]};

%%% Plot stable manifold (estimated via a straight line from the eigenvectors of Jacobian near saddle point)
x1 = 0:0.01:0.6;
figure(1),
plot(x1,2.12194-3.44874*x1,'k--','LineWidth',2),hold on

%%% Plot steady state points with filled circles signifying stable s.s and triangle as saddle point
figure(1),
plot(0.27189,2.447,'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(0.054,0.486,'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(0.17045,1.5341,'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on

for i=1:length(hist)
    C = {['r','b']}; % Assign colors for trajectory plot

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol = dde23(@ddefunc, tau, hist{i}, [Tstart Tend], opts);

    %%% Solving ODE
    solODE = ode23s(@odefunc, [Tstart Tend], hist{i}, []);

    %%% Plot trajectory
    figure(1),
    plot(sol.y(1,:),sol.y(2,:),'LineWidth',1,'Color',C{1}(i)), hold on
    plot(solODE.y(1,:),solODE.y(2,:),'k','LineWidth',2), hold on
    xlabel('G')
    ylabel('P')
    box off
    set(gca,'FontSize',20)
end
figure(1), axis([0 0.4 0 2.6])
xticks([0 0.2 0.4])
yticks([0 1.3 2.6])

%% Figure 8C--------------------------------------------------------------

%%% Solving DDE
opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
sol2 = dde23(@ddefunc, tau, hist{1}, [Tstart Tend],opts);

%%% Plotting timeseries
figure(2),
plot(sol2.x/(2*tau),sol2.y(2,:),'Color',[0 0.6 0],'LineWidth',2)
ylabel('P')
xlabel('oscillation number')
box off
set(gca,'FontSize',20)
axis([0 60 0 3])

%------------------------------------------------------------------------

function dout = ddefunc(t, y, yl)
%%% Differential equations function

% Define the variables
g = y(1);
p = y(2);
glag = yl(1,:);
plag = yl(2,:);

% Parameters that give bistability
alpha = 0.1;
gamma = 2;
k = 2;
n = 3;
lambda = 9;
delta = 1;

% Differential Equation
dg = alpha*(1-g)-gamma*g+(1-g)*(plag^n/(k^n+plag^n));
dp = lambda*glag-delta*p;

dout = [dg dp]';
end

%------------------------------------------------------------------------

function dout = odefunc(t, y)
%%% Differential equations function

% Define the variables
g = y(1);
p = y(2);

% Parameters 
alpha = 0.1;
gamma = 2;
k = 2;
n = 3;
lambda = 9;
delta = 1;

% Differential Equation
dg = alpha*(1-g)-gamma*g+(1-g)*(p^n/(k^n+p^n));
dp = lambda*g-delta*p;

dout = [dg dp]';
end
