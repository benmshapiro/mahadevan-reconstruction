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
v = [0;1;0;0];         % initial flux [ac+o2, gluc+o2, glu+o2->ac, glu->ac]
                       % assuming only glu+o2 is active at the begining
X = 0.001;             % initial biomass in ??? unit

% Define constants
A = msm.S;             % 4x4 stoichiometric matrix of simplfied network (fig 2)
lb = zeros(4,1);       % lower bounds of 4 rxns of ac+o2, gluc+o2, glu+o2->ac, glu->ac
kla = 7.5;  %O2 mass transfer coefficient
our = 15;   %O2 uptake constraint
ogp = 0.21; %O2 gas phase concentration 
Km = 0.015; %Glucose saturation constant
vdot = [0.1; 0.3; 0.3; 0.1]; % max rates of change in flux rate per time step for 4 rxnx
weights = -1*[1; 1; 1; 1];   % weights of 4 rxns in optimization
deltaT = 1/1000;             % time step in 1/hr
length = 10000;              % number of time step
metaProf = zeros(length,4);  % concentrations of metabolites (glc, ac, o2, X) 
metaProf(1,:) = [z' X];
fluxes = zeros(length,4);
fluxes(1,:) = zeros(1,4);

%% Model Simulation

for ii = 1:length
   vold = v;             % 4x1 vector
   zold = z;             % 3x1 vector
   Xold = X;             % 1x1
   
   % Formulate and Solve LP
   AnonNegMets = -1*A*deltaT;
   bnonNegMets = [zold; Xold];               % initial conditions (metabolites & biomass),vector 4x1
   
   AlimVchange1 = eye(4);                    % 4x4 matrix
   blimVchange1 = (vold + vdot);             % max flux of 4 rxns, vectoro 4x1
   
   AlimVchange2 = -1*eye(4);                 % identity maxtrix
   blimVchange2 = (-1*vold + vdot);          % ??? minimum uptake rate, vector 4x1
   
   AglucUptake = -1*A(1,:);                  % stoichiometric coeff of glu in 4 rxns, 4x1 vector
   bglucUptake = 10*zold(1)/(Km + zold(1));
   Ao2uptake = -1*A(3,:);                    % 4x1 vectore
   bo2uptake = our;                          % o2 uptake is fixed
      
   if ii <= 1
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
interval = 0.001:0.001:10;

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