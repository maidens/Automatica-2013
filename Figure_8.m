%% Approximation of the viability kernel using support functions and support vectors
%     for pharmacokinetic example
%
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires CVX, Level Set Toolbox, Ellipsoidal Toolbox, Multi-Parametric
%     Toolbox


clear all; clc
cvx_quiet(true);

% note: times in minutes, drug masses in mcg, volumes in mL


%% Paedfusor model parameters (ages 1-12 years)

weight = 35;    % patient weight in kg

k10=0.1527*weight^-0.3;
k12=0.114;
k13=0.0419;
k21=0.055;
k31=0.0033;
keo=0.26;       % not used - has to do with effect site concentration
V1=458.4*weight;
V2=V1*k12/k21;  % not used 
V3=V1*k13/k31;  % not used


%% continuous-time state matrices
A = [-(k10+k12+k13)  k12  k13;
           k21      -k21   0 ;
           k31        0  -k31 ];
  
B = [1/V1;
      0  ;
      0  ];

C = [1 0 0];

D = 0;


%% continuous-time pade approximation

sys = ss(A,B,C,D,'InputDelay',30/60); % 30 second delay
sys_pade = pade(sys,3);

A_pade = sys_pade.A
B_pade = sys_pade.B

n = length(sys_pade.A); % state space dimension


%% discretization
tau = 90;            % time horizon (minutes)
steps = 4*90;        % number of time steps for support function method
rho = tau/steps;     % disctetization interval for support function method

sysd = c2d(sys_pade,rho);

Ad = sysd.A;
Bd = sysd.B;
Astar = inv(Ad)';


%% safe set is image of unit ball under x |-> (A_K0 x + q_K0)
%    units of mcg/mL

A_K0 = diag([2.5 5 5 100 100 100]); 
q_K0 = [3.5; 5; 5 ; 0; 0; 0];


%% input constraints are [0, u_max] (in u-space)
%    (i.e. image of unit ball under x |-> (A_V x + q_V) in v-space)

max_input = 200;           % in mcg/kg/min
u_max = weight*max_input;  % converted to mcg/min        

A_V = -u_max/2*Bd;
q_V = -u_max/2*Bd;


%% Setup approximation directions 

% vertices of an icosahedron 
phi=(1+sqrt(5))/2;
V1=[0;0;0;0;-1;-1;1;1;-phi;phi;phi;-phi;];
V2=[-1;-1;1;1;-phi;phi;phi;-phi;0;0;0;0;];
V3=[-phi;phi;phi;-phi;0;0;0;0;-1;-1;1;1;];
V = [V1 V2 V3];      % vertices of an icosahedron 

L = [      V'        zeros(3) zeros(3)
     zeros(size(V'))  eye(3)  -eye(3) ];
L = A_K0\L;          % approximation directions
num_L = size(L,2);   % number of support vectors (approximation directions)


%% Run support vector method to compute over-approximating polytope
%   (A_viab, b_viab) and under-approxiamting polytope v for viability kernel

tic

computeViab_SupportVect

time_8 = toc


%% construct geometric objects using ellipsoidal and multi-parametric
%    toolboxes

v_proj = v(:,1:3);
U = A_K0(1:3,1:3)*ell_unitball(3)+q_K0(1:3); % ellipsoidal viable set (using ET)
P = polytope(A_viab,b_viab);
P_over_7000 = projection(P,[1 2 3]);   % polytopic over-approxiamtion (MPT)
P_under_7000 = polytope(v_proj);       % polytopic under-approximation (MPT)


%% plot the above objects

Options1.color=2*[0 0.1569 0.3490];        % UBC blue
Options2.color=[0.4549 0.5686 0.6392];     % UBC grey
OptionsU1.shade=0.25;
Options1.shade=1;
Options2.shade=0.3;
OptionsU2.shade=1;

figure                     
hold on

    plot(U,OptionsU1)
    plot(P_under_7000,Options1)
 
    light('Position',[10 0 10],'Style','infinite');
    light('Position',[0 0 10],'Style','infinite');
    light('Position',[0 12 6],'Style','infinite');
    light('Position',[0 -6 0],'Style','infinite');
    
hold off
axis([0 6 0 10 0 10])

figure                     
hold on

    plot(U,OptionsU2)
    plot(P_over_7000,Options2)
 
    light('Position',[10 0 10],'Style','infinite');
    light('Position',[0 0 10],'Style','infinite');
    light('Position',[0 12 6],'Style','infinite');
    light('Position',[0 -6 0],'Style','infinite');
    
hold off
axis([0 6 0 10 0 10])





    
    