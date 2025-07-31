function [cost,sol]=eadymodel(K,c,eadyparams)

% Integrates the Eady model in Zurita-Gotor and Held (2025) and evaluates the cost function. 
%
% [cost,sol]=eadymodel(K,c,[eadyparams])
%
% Input arguments are nondimensional wavenumber (a real) and phase speed (a complex, in general)
% 
% The optional input argument eadyparams should be of the form {L, Lambda} where L is the jet latitude 
% and Lambda the poleward shear. When either of these are not indicated, defaults are taken from the
% control case in the paper {3.5, -0.5}
%
% Output arguments are the value of the cost function (cost, defined by Eq. 15 in the paper) and  
% the full solution to the ODE (sol). This is a structure consisting of the following fields:
%     sol.y: latitude (a list of real values)
%     sol.u: zonal velocity (a list of complex values)
%     sol.v: meridional velocity (a list of complex values)
%     sol.z: geopotential (a list of complex values)

%Read model parameters, using defaults when necessary
if ~exist('eadyparams') eadyparams={}; end  %Model parameters not input
if length(eadyparams)<2  Lambda=-0.5;       %Default poleward shear
else Lambda=eadyparams{2}; end
if length(eadyparams)<1  L=3.5;             %Default jet latitude 
else L=eadyparams{1};end

%Define and integrate ODE (Eq. 9 in the paper)
udot=@(y,u) [u(2); K^2*(1-(.5*y^2-c)^2)*u(1)];
[y, u]=ode45(udot,[0 L],[1 0]);

%Evaluate cost function
U=.5*L^2;                  %Nondimensional velocity at the front
m_over_L=sqrt(1+K^2/L^2);  %Use 1 for longwave approximation
cost=(U-c)*(1+m_over_L*(U-c))*u(end,1)-u(end,2)*(L-Lambda)/K^2;  %Eq. 15

if nargout>1 %For efficiency when iterating, calculate full solution only when requested 
  sol.y=y;                     %latitudinal coordinate
  sol.u=u(:,1);                %zonal velocity
  sol.v=-sqrt(-1)*u(:,2)/K;    %meridional velocity
  sol.z=-(.5*y.^2-c).*u(:,1);  %geopotential
end        


end


