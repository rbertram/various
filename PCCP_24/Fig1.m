% Fig1.m

% Used to make Figure 1B-C and E-F from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% This code solves simple ODE systems modelling both mutual negative feedback
% and mutual positive feedback of 2 proteins. After, phase-plane and time series are plotted for various initial conditions.

% Parameters_NegNeg.m, Parameters_PosPos, SteadyState_Solutions_NegNeg.m,
% and SteadyState_Solutions_PosPos.m must be in current directory!!

% X -| Y inhibition of Y by X
% Y -| X inhibition of X by Y
% ODE Equations:
% dx/dt = alphaX*y^n / (1 + y^n) - x;
% dy/dt = alphaY*x^n / (1 + x^n) - y;

% X -| Y inhibition of Y by X
% Y -| X inhibition of X by Y
% ODE Equations:
% dx/dt = alphaX / (1 + y^n) - x;
% dy/dt = alphaY / (1 + x^n) - y;

% Author: Christopher Ryzowicz
% -----------------------------------------------------------------------
clear all
clc
close all
% -------------------------------------------------------------------------
%%% Time interval
Tstart = 0;
Tend = 15;

%%% Parameters Neg-Neg
Parameters_NegNeg

%%% Parameters Pos-Pos
Parameters_PosPos

%%% Steady State Solutions Neg-Neg
[ss_x, ss_y] = SteadyState_Solutions_NegNeg();
ss_x = real(ss_x(imag(ss_x) == 0)); % Extract only the real steady states
ss_y = real(ss_y(imag(ss_y) == 0)); % Extract only the real steady states
%%% Steady State Solutions Pos-Pos
[ss_x2, ss_y2] = SteadyState_Solutions_PosPos();
ss_x2 = real(ss_x2(imag(ss_x2) == 0)); % Extract only the real steady states
ss_y2 = real(ss_y2(imag(ss_y2) == 0)); % Extract only the real steady states

%%% Generating points for stable manifold of Neg-Neg
x1 = 0:0.01:10;
%%% Generating points for stable manifold of Pos-Pos
x2 = 0:0.01:2.15;

%%% Plot the stable manifold for each system
figure(1),
plot(x1,x1,'k--','LineWidth',2),hold on
figure(2),
plot(x2,1.36-x2,'k--','LineWidth',2),hold on

%%% Plot steady state points with filled circles signifying stable s.s and triangle as saddle point
figure(1),
plot(ss_x(2:3),ss_y(2:3),'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
plot(ss_x(1),ss_y(1),'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'),hold on
figure(2),
plot(ss_x2(1:2:3),ss_y2(1:2:3),'o','Color','k','MarkerSize',16,'MarkerFaceColor','k'), hold on
plot(ss_x2(2),ss_y2(2),'^','Color','k','MarkerSize',16,'MarkerFaceColor','k'), hold on

%%% History constants (x_0,y_0)
histNeg = {[9,7], [4,5], [8,3], [2,6], [1,0.8], [1,1.5]};
histPos = {[0.25,0.75], [1.0,0.6], [1.75,1], [0.75,0.2], [0.5,1.6], [0.45,0.85]};

for i=1:length(histNeg)

    C = {['r','g','b','c','m','k']}; % Specify color order for plots

    %%% Solve ODE systems
    solNeg = ode23s(@odefunc1, [Tstart Tend], histNeg{i}, [], p);
    solPos = ode23s(@odefunc2, [Tstart Tend], histPos{i}, [], q);

    %%% Plot Neg-Neg phase plane
    figure(1),
    plot(solNeg.y(1,:),solNeg.y(2,:),'LineWidth',2,'Color',C{1}(i)), hold on
    xlabel('x')
    ylabel('y')
    title('Phase plane for mutual inhibition model')
    box off
    set(gca,'FontSize',20)
    %%% Plot Pos-Pos phase plane
    figure(2),
    plot(solPos.y(1,:),solPos.y(2,:),'LineWidth',2,'Color',C{1}(i)), hold on
    xlabel('x')
    ylabel('y')
    title('Phase plane for mutual activation model')
    box off
    set(gca,'FontSize',20)
    %%% Plot Neg-Neg time series
    figure(3),
    subplot(2,1,1),plot(solNeg.x,solNeg.y(1,:),'Color',C{1}(i),'LineWidth',2), hold on
    ylabel('x')
    title('Time series for mutual inhibition model')
    box off
    set(gca,'FontSize',20)
    subplot(2,1,2),plot(solNeg.x,solNeg.y(2,:),'Color',C{1}(i),'LineWidth',2), hold on
    xlabel('t')
    ylabel('y')
    box off
    set(gca,'FontSize',20)
    %%% Plot Pos-Pos time series
    figure(4),
    subplot(2,1,1),plot(solPos.x,solPos.y(1,:),'Color',C{1}(i),'LineWidth',2), hold on
    ylabel('x')
    title('Time series for mutual activation model')
    box off
    set(gca,'FontSize',20)
    subplot(2,1,2),plot(solPos.x,solPos.y(2,:),'Color',C{1}(i),'LineWidth',2), hold on
    xlabel('t')
    ylabel('y')
    box off
    set(gca,'FontSize',20)

end
figure(1), axis([0 10.5 0 10.5])
xticks([0 5 10])
yticks([0 5 10])
figure(2), axis([0 2 0 2])
xticks([0 1 2])
yticks([0 1 2])
for i =1:2
    figure(3),subplot(2,1,i),
    axis([0 15 0 10])
    xticks([0 5 10 15])
    yticks([0 5 10])
    figure(4),subplot(2,1,i),
    axis([0 15 0 2])
    xticks([0 5 10 15])
    yticks([0 1 2])
end
%------------------------------------------------------------------

function doutNeg = odefunc1(t, out, p)
%%% Define the variables
x = out(1);
y = out(2);

%%% Differential Equation
dx = p.alphaX / (1 + y^p.n) - x;
dy = p.alphaY / (1 + x^p.n) - y;

doutNeg = [dx dy]';
end

function doutPos = odefunc2(t, out, q)
%%% Define the variables
x = out(1);
y = out(2);

%%% Differential Equation
dx = (q.alphaX*y^q.n) / (1 + y^q.n) - x;
dy = (q.alphaY*x^q.n) / (1 + x^q.n) - y;

doutPos = [dx dy]';
end
