% SteadyState_Solutions_PosPos.m

function [solx, soly] = SteadyState_Solutions_PosPos()

syms x y
% import parameters from Parameters_PosPos.m file
Parameters_PosPos

[solx,soly] = solve((q.alphaX*y^q.n)/(1+y^q.n)-x==0, (q.alphaY*x^q.n)/(1+x^q.n)-y==0);

end
