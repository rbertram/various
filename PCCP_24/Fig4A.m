% Fig4A.m

% Used to make Figure 4A from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% This code simulates the negative-negative DDE and calculates how long the oscillations
% last as a function of the delay, tau. The method of calculating how long these oscillations last is given in the text of the
% manuscript.

% Parameters_NegNeg.m and SteadyState_Solutions_NegNeg.m must be in current directory!!

% System of equations for negative-negative DDE model:
% dx/dt = alphaX / (1 + ylag^n) - x;
% dy/dt = alphaY / (1 + xlag^n) - y;

% Author: Christopher Ryzowicz 
%-------------------------------------------------------------------------
clear all
clc
%  clf
%-------------------------------------------------------------------------

% Solve for steady state values for negative-negative model
[ss_x1, ss_y1] = SteadyState_Solutions_NegNeg();
% Caculate midpoint of bistable steady state values
ss_midpoint1 = double((max(ss_x1)+min(ss_x1))/2);

%-------------------------------------------------------------------------

%%% Initialize list of delay values
tau_list = 2:0.5:8; %10

for i=1:length(tau_list)
i
    %%% Timespan and delay
    tau = tau_list(i); % Delay values
    Tstart = 0; % Start time of solution
    Tend = 40000; % End time of solution

    %%% Parameters
    Parameters_NegNeg

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol1 = dde23(@ddefunc1, tau, @yhist, [Tstart Tend],opts,p); % Solution for negative-negative model

    %%% Assign beta weight value for convex combination depending on if x reaches high or low steady state
    %%% Calculate time required for last 'x' value of oscillation to reach convex combination value with weight beta
    if sol1.y(1,end) < ss_midpoint1
        beta = 0.1;
        convex_combo = (1-beta)*min(ss_x1) + beta*max(ss_x1);
        oscillate_timeX1(i) = sol1.x(find(sol1.y(1,:)>convex_combo,1,'last'));
    else
        beta = 0.9;
        convex_combo = (1-beta)*min(ss_x1) + beta*max(ss_x1);
        oscillate_timeX1(i) = sol1.x(find(sol1.y(1,:)<convex_combo,1,'last'));
    end
end

%%% Save data to a .mat file
% save('CodeXX8_Data.mat')

%------------------------------------------------------------------------

%%% Plot scatter of oscillation duration vs. tau value
figure(2),scatter(tau_list,oscillate_timeX1./(2*tau_list),'filled','SizeData',100)
xlabel('\tau')
ylabel('Number of Oscillations')
box off
set(gca,'FontSize',20)

%%% Save data to a .mat file
% save('CodeXX8_Data.mat')

%----------------------------------------------------------------------

function dout = ddefunc1(t, out, yl, p)

%%% Define the variables
x = out(1);
y = out(2);
xlag = yl(1,:);
ylag = yl(2,:);

%%% Differential Equation
dx = p.alphaX / (1 + ylag^p.n) - x;
dy = p.alphaY / (1 + xlag^p.n) - y;

dout = [dx dy]';
end

%------------------------------------------------------------------------

function y = yhist(t, p)
%%% History function

x_0 = 4;
y_0 = 3;

y = [x_0 y_0];
end
