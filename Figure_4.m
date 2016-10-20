%% Exact computation of the viability kernel using polytopes
%
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires Multi-Parametric Toolbox

clear all; clc

%% Continuous-time double integrator

A = [0 1;
     0 0];
 
B =  [ 0 ;
       1  ]; 
   
C = [1 0];

D = 0;
   
sys = ss(A,B,C,D);

U = unitbox(1,0.3);   % Input constraint set

K0 = unitbox(2,0.5);  % State constraint set


%% Discretization

N = 40;     % maximum number of time steps
tau = 4/N;  % time step (2 second horizon)

sys_mpt = mpt_sys(sys,tau);

% Inversion of time direction (to get backward reach sets)
sys_mpt.A = inv(sys_mpt.A); 
sys_mpt.B = -inv(sys_mpt.A)*sys_mpt.B;

% Required arguments for MPT but not used for reach computations
sys_mpt.xmax = [.5 .5]'; sys_mpt.xmin = [-.5 -.5]'; 
sys_mpt.umax = .3; sys_mpt.umin = -.3;



%% Run Algorithm 3

tic
computeViab_Polytope
time = toc


%% Plot the results

Options1.color=[0 0.1569 0.3490];      % UBC blue
Options2.color=[0.4549 0.5686 0.6392]; % UBC grey

figure
    hold on
    plot(K0,Options2)
    plot(KN,Options1)
    hold off
    axis([-.6,.6,-.6,.6])
    set(gca,'DataAspectRatio',[1 1 1])


