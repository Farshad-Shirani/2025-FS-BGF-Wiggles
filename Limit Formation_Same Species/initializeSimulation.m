function [initialPopulations, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation

%---simulation parameters----------------------------------------------
T = 1500; %final time
storageTimes = 0 : 2 : T; % solution storage times. Solutins are sored at closest time samples depending on Dt
x_0 = 0;    x_I = 50; % x1 range
simulationParameters = struct('T', T, 'times', storageTimes, 'x_0', x_0, 'x_I', x_I);

%---discretization parameters------------------------------------------
Dt = 0.005;  % time steps
Dx = (x_I - x_0)/ 500; %x1 mesh   
discretizationParamaters = struct('Dt', Dt, 'Dx', Dx);

x = x_0 : Dx : x_I;
I = length(x);

%---model parameters---------------------------------------------------
global epsilon 
global AlleeFlag J1 J2
global D1 D2 V_s V_u1 V_u2 U kappa K1 K2 Q_opt R1 R2

epsilon = 1e-5; %d/dx(log (N + epsilon)) % (to avoid singularity due to the d/dx(log N) term

AlleeFlag = 1; % 1: with Allee effect, 0: without Allee effect
J1 = 0.02 * ones(I, 1);  % critical density for population 1
J2 = 0.02 * ones(I, 1);  % critical density for population 2

D1 = 1;  % diffusion matrix for populationn 1
D2 = 1;  % diffusion matrix for populationn 2
V_s = 1/0.2;    % 1/measure of the strength of stabilizing selection (V_s = 1/S, where S defined in Table 1)
V_u1 = 9;    % variance of the within phenotipic-resource utility curve for population 1 (V_u1 = V_1 where V is defined in Table 1)
V_u2 = 9;    % variance of the within phenotipic-resource utility curve for population 2 (V_u2 = V_2 where V is defined in Table 1)
U = 0.02;  % rate of increase in trait variance due to mutation
kappa = 0;  % population impact factor
K1 = 1 * ones(I, 1);  % carrying capacity for population 1
K2 = 1 * ones(I, 1);  % carrying capacity for population 2
R1 = 1 * ones(I, 1);  % intrinsic rate of increase for optimally adapted individuals of population 1
R2 = 1 * ones(I, 1);  % intrinsic rate of increase for optimally adapted individuals of population 2

dQ_opt = -0.5; % optimum trait gradient dQ
Q_0 = 10; % optimum trait at x_J

Q_opt = Q_0 + dQ_opt  * (x - x(end))';    % optimum trait value 

%-----------------

modelParameters = struct('epsilon', epsilon, 'AlleeFlag', AlleeFlag, 'J1', J1, 'J2', J2, ...
    'D1', D1, 'D2', D2,  'V_s', V_s, 'V_u1', V_u1, 'V_u2', V_u2, 'U', U, 'kappa', kappa, ...
    'K1', K1, 'K2', K2, 'Q_opt', Q_opt, 'R1', R1, 'R2', R2);

%---solver parameters----------------------------------------------------
eta = 1/2;

%---initial values-------------------------------------------------------
initialPopulations = struct('density', [], 'trait_mean', [], 'trait_variance', []);

%---For same species --------------------------------------------
center = 15; % population center
radius = 2; % population radius in each direction 
%initialPopulations(1).density = 0.5 * sech( abs( x - center ) / radius )';
initialPopulations(1).density = 0.5 * bump1( ( x - center ) / radius )';
center_index_x = dsearchn(x', center); 
initialPopulations(1).trait_mean = 0.1 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
initialPopulations(1).trait_variance = 1 * ones(I, 1);

center = 35; % population center 
radius = 2; % population radius in each direction
%initialPopulations(2).density = 0.5 * sech( abs( x - center ) / radius )';
initialPopulations(2).density = 0.5 * bump1( ( x - center ) / radius )';
center_index_x = dsearchn(x', center); 
initialPopulations(2).trait_mean = 0.1 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
initialPopulations(2).trait_variance = 1 * ones(I, 1);

end

