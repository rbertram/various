% Fig3A.m

% Used to make Figure 3B-C from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.


% This code solves a simple DDE system modelling mutual negative feedback of 2
% molecules. After, we plot 2 trajectories, one from ODE and other from
% DDE. In addition, we provide a sample time series in x for the DDE system.

% Parameters_NegNeg.m and SteadyState_Solutions_NegNeg.m must be in current directory!!

% X -| Y inhibition of Y by X
% Y -| X inhibition of X by Y
% Delay Equation:
% dx/dt = alphaX / (1 + ylag^n) - x;
% dy/dt = alphaY / (1 + xlag^n) - y;

% Author: Christopher Ryzowicz 
%-------------------------------------------------------------------------
clear all
clc
close all
%% Figure 3B--------------------------------------------------------------
%%% Timespan and delay
tau = 6;% Delay value
Tstart = 0; % Start time of solution
Tend = 600; % End time of solution

%%% Parameters
Parameters_NegNeg

%%% Steady State Solutions
[ss_x, ss_y] = SteadyState_Solutions_NegNeg();
ss_x = real(ss_x(imag(ss_x) == 0)); % Extract only the real steady states
ss_y = real(ss_y(imag(ss_y) == 0)); % Extract only the real steady states

%%% History constants (x_0,y_0)
hist = {[4,3], [3.5,5]};

%%% Plot the stable manifold (stable manifold for this symmetric system is the 'y=x' line)
x1 = 0:0.01:10;
figure(1),
plot(x1,x1,'k--','LineWidth',2),hold on

%%% Plot steady state points with filled circles signifying stable s.s and triangle as saddle point
figure(1),
plot(ss_x(2:3),ss_y(2:3),'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'), hold on
plot(ss_x(1),ss_y(1),'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'), hold on

for i=1:length(hist)
    C = {['r','b']}; % Assign colors to trajectory plot

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol = dde23(@ddefunc, tau, hist{i}, [Tstart Tend], opts, p);

    %%% Solving ODE
    solODE = ode23s(@odefunc, [Tstart Tend], hist{i}, [], p);

    %%% Plot trajectory of ODE soln
    figure(1),
    plot(sol.y(1,:),sol.y(2,:),'LineWidth',1,'Color',C{1}(i)), hold on
    plot(solODE.y(1,:),solODE.y(2,:),'k','LineWidth',2), hold on
    xlabel('x')
    ylabel('y')
    box off
    set(gca,'FontSize',20)
end

%%% Adjust axes
figure(1), axis([0 10.5 0 10.5])
xticks([0 5 10])
yticks([0 5 10])

%% Figure 3C-----------------------------------------------------------

%%% Solving DDE
opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
sol2 = dde23(@ddefunc, tau, hist{1}, [Tstart Tend],opts,p);

%------------------------------------------------------------------------

%%% Plot individual state variables against oscillation number (divide the time by period of oscilllation which here is approximately 2*tau)
figure(2),
plot(sol2.x/(2*tau),sol2.y(1,:),'Color',[0 0.45 0.73],'LineWidth',2), hold on
xlabel('oscillation number')
ylabel('x')
yticks([0 2 4 6 8 10])
box off
set(gca,'FontSize',20)

%----------------------------------------------------------------------

function dout = ddefunc(t, out, yl, p)
%%% Define the variables
x = out(1);
y = out(2);
xlag = yl(1);
ylag = yl(2);

%%% Differential Equation
dx = p.alphaX / (1 + ylag^p.n) - x;
dy = p.alphaY / (1 + xlag^p.n) - y;

dout = [dx dy]';
end

function dout1 = odefunc(t, out, p)
%%% Define the variables
x = out(1);
y = out(2);

%%% Differential Equation
dx = p.alphaX / (1 + y^p.n) - x;
dy = p.alphaY / (1 + x^p.n) - y;

dout1 = [dx dy]';
end
