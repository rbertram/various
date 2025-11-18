% Fig5B.m

% Used to make Figure 5B from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% This code simulates the positive-positive DDE model and calculates how long the oscillations
% last as a function of distance from the stable manifold (seperatrix). The method of calculating how long these oscillations last is given in the text of the
% manuscript.

% Parameters_PosPos and SteadyState_Solutions_PosPos.m must be in current directory!!

% System of equations for positive-positive DDE model:
% dx/dt = (alphaX*ylag^n) / (1 + ylag^n) - x;
% dy/dt = (alphaY*xlag^n) / (1 + xlag^n) - y;

% Author: Christopher Ryzowicz
%-------------------------------------------------------------------------
clear all
clc
%  clf
%-------------------------------------------------------------------------

% Solve for steady state values for positive-positive model
[ss_x2, ss_y2] = SteadyState_Solutions_PosPos();
% Caculate midpoint of bistable steady state values
ss_midpoint2 = double((max(ss_x2)+min(ss_x2))/2);

%-------------------------------------------------------------------------

% Initialize history function pairings
x0_list = 0.82:0.005:0.92; % This vector gives points above stable manifold
% x0_list = 0.78:-0.005:0.68; % This vector gives points below the stable manifold
y0_list = x0_list - 0.24;
hist_list = {x0_list, y0_list};

% Calculate distance list to stable manifold using predetermined orthoganal line y=x-0.24 and
% history values in reference to point (0.8,0.56)
dist_to_manifold = sqrt((x0_list-0.8).^2 + (y0_list-0.56).^2);

for i=1:length(hist_list{1})
i
    %%% Timespan and delay
    tau = 6; % Delay values
    Tstart = 0; % Start time of solution
    Tend = 2000; % End time of solution

    %%% Parameters
    Parameters_PosPos

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol2 = dde23(@ddefunc2, tau, [hist_list{1}(i), hist_list{2}(i)], [Tstart Tend],opts,q); % Solution for positive-positive model

    %%% Assign beta weight value for convex combination depending on if x reaches high or low steady state
    %%% Calculate time required for last 'x' value of oscillation to reach
    %%% convex combination value with weight beta
    if sol2.y(1,end) < ss_midpoint2
        beta = 0.1;
        convex_combo = (1-beta)*min(ss_x2) + beta*max(ss_x2);
        oscillate_timeX2(i) = sol2.x(find(sol2.y(1,:)>convex_combo,1,'last'));
    else
        beta = 0.9;
        convex_combo = (1-beta)*min(ss_x2) + beta*max(ss_x2);
        oscillate_timeX2(i) = sol2.x(find(sol2.y(1,:)<convex_combo,1,'last'));
    end
end

%------------------------------------------------------------------------

%%% Plot scatter of oscillation duration vs. distance to stable manifold
figure(2),scatter(dist_to_manifold,oscillate_timeX2./(2*tau),'filled','SizeData',100)
xlabel('\Delta')
ylabel('Number of Oscillations')
box off
set(gca,'FontSize',20)

%%% Save data as a .mat file
% save('CodeXX6_Data.mat')

%----------------------------------------------------------------------

function dout = ddefunc2(t, out, yl, q)

%%% Define the variables
x = out(1);
y = out(2);
xlag = yl(1,:);
ylag = yl(2,:);

%%% Differential Equation
dx = (q.alphaX*ylag^q.n) / (1 + ylag^q.n) - x;
dy = (q.alphaY*xlag^q.n) / (1 + xlag^q.n) - y;

dout = [dx dy]';
end
