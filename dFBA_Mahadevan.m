% Description: Replicating Mahadevan et al. (2002) Bioph. Journ. --> SOA formalism
% Author: Benjamin Shapiro, University of Oregon
% Date: 2015-01-19
% Modified: 2015-02-12
% Comments: Use in conjunction with MahadevanSimpleModel.m
%
%% Model Steup
% Load simple metabolic model
MahadevanSimplifiedModel;

% Initialize variables
z = [10.8; 0.4; 0.21]; % Glucose, Acetate, Oxygen (all in mM)
v = [0;1;0;0];
X = 0.001;

% Define constants
A = msm.S;
lb = zeros(4,1);
kla = 7.5;  %O2 mass transfer coefficient
our = 15;   %O2 uptake constraint
ogp = 0.21; %O2 gas phase concentration 
Km = 0.015; %Glucose saturation constant
vdot = [0.1; 0.3; 0.3; 0.1]; %flux rate of change constraints
weights = -1*[1; 1; 1; 1];
deltaT = 1/1000;
length = 5000;
metaProf = zeros(length,4); % metabolite profiles
metaProf(1,:) = [z' X];
fluxes = zeros(length,4);
fluxes(1,:) = zeros(1,4);

%% Model Simulation

for ii = 1:length
   vold = v;
   zold = z;
   Xold = X;
   
   % Formulate and Solve LP
   AnonNegMets = -1*A*deltaT;
   bnonNegMets = [zold; Xold];
   AlimVchange1 = eye(4);
   blimVchange1 = (vold + vdot);
   AlimVchange2 = -1*eye(4);
   blimVchange2 = (-1*vold + vdot);
   AglucUptake = -1*A(1,:);
   bglucUptake = 10*zold(1)/(Km + zold(1));
   Ao2uptake = -1*A(3,:);
   bo2uptake = our;
      
   if ii < 1
        Aineq = [ AnonNegMets; AlimVchange1; AlimVchange2; AglucUptake; Ao2uptake];
        bineq = [ bnonNegMets; blimVchange1; blimVchange2; bglucUptake; bo2uptake];
   else
        Aineq = [ AnonNegMets; AglucUptake; Ao2uptake];
        bineq = [ bnonNegMets; bglucUptake; bo2uptake];
   end
   
   v = cplexlp(weights, Aineq,bineq,[],[],lb,[]);

   % Integrate to find new extracellular concentrations
   dGluc = A(1,:)*v*Xold;
   dAcet = A(2,:)*v*Xold;
   dO2 = A(3,:)*v*Xold + kla*(0.21-zold(3));
   dX = sum(v)*X;
   
   z = [dGluc; dAcet; dO2]*deltaT + zold;
   X = dX*deltaT + Xold;
   metaProf(ii,:) = [z' X];
   fluxes(ii,:) = v';
   
end
%% Post-Processing

% Plotting
interval = 0.001:0.001:5;

subplot(2,2,1);plot(interval, fluxes(:,1),interval, fluxes(:,2),interval, fluxes(:,3),interval, fluxes(:,4));
title('Flux rate');
legend('v1','v2','v3','v4');
xlabel('Time(hr)');
ylabel('mmol g^-^1 hr^-^1');

subplot(2,2,2);[AX,H1,H2] = plotyy(interval, metaProf(:,2),interval, metaProf(:,4));
title('Acetate and Biomass');
set(AX(1),'YLim',[0 10]);
set(AX(1),'YTick',0:2:10)
set(get(AX(1),'Ylabel'),'String','Acetate conc. (mM)');
set(get(AX(2),'Ylabel'),'String','Biomass rate (hr^-^1)'); 
xlabel('Time(hr)');

subplot(2,2,3);plot(interval, metaProf(:,3));
title('Oxygen');
ylabel('Oxygen conc.(mM)');
xlabel('Time(hr)');

subplot(2,2,4);plot(interval, metaProf(:,1));
title('Glucose');
ylabel('Glucose conc.(mM)');
xlabel('Time(hr)');