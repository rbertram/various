% Fig7.m

% Used to make Figure 7A-B from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% This code solves an assymetric version of delayed mutual inhibiton or
% delayed mutual activation between 2 species depending on the sign of the
% cooperativity parameters n1,n2. Time series of 'x' species is plotted.
% Feel free to change delay values (line 27), cooperativity parameter (line 33), or initial
% history functions (lines 67-74)to explore the dynamics.

% Figure 7A
% X -| Y inhibition of Y from X after delay tau_1
% X |- Y inhibition of X from Y after delay tau_2
% Figure 7B
% X -> Y activation of Y from X after delay tau_1
% X <- Y activation of X from Y after delay tau_2

% System of equations:
% dx/dt = alpha / (1 + y(t-tau1)^n1) - x;
% dy/dt = alpha / (1 + x(t-tau2)^n2) - y;
% n1,n2 > 0: mutual inhibition
% n1,n2 < 0: mutual activation

%   Author: Christopher Ryzowicz

%-------------------------------------------------------------------------

%%% Timespan and delay
tau = [8 12];% Delay values [tau1, tau2]
Tstart = 0; % Start time
Tend = 2000; % End time 

%%% Parameters
p.alpha = 10;
p.n = [2 2]; %[n1, n2]

%%% Solving DDE
opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
sol = dde23(@ddefunc, tau, @yhist, [Tstart Tend],opts,p);

%------------------------------------------------------------------------


%%% Plot individual timeseries
figure(2), 
plot(sol.x,sol.y(1,:))
xlabel('time')
ylabel('x')
box off
set(gca,'FontSize',20)

%----------------------------------------------------------------------
function dout = ddefunc(t, out, yl, p)
% Define the variables
x = out(1);
y = out(2);
xlag = yl(1,:);
ylag = yl(2,:);

% Differential Equation
dx = p.alpha / (1 + ylag(1)^p.n(1)) - x;
dy = p.alpha / (1 + xlag(2)^p.n(2)) - y;

dout = [dx dy]';
end

%------------------------------------------------------------------------

function y = yhist(t, p)
%%% History function

x_0 = 7;
y_0 = 3;

y = [x_0 y_0];
end
