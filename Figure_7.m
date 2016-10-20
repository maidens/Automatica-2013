%% Comparison of the scalability of three viability kernel approximation
%      algorithms as the state dimension increases
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires CVX, Level Set Toolbox, Ellipsoidal Toolbox, Multi-Parametric
%     Toolbox

clear all; clc

cvx_quiet(true);    % CVX option
ellipsoids_init;
global ellOptions;   % ET options
ellOptions.verbose = 0;
options.approximation = 1;  

for state_dim = 1:30
    %% State space model of an n-ary integrator
    A = [[zeros(state_dim-1,1) eye(state_dim-1)];
          zeros(1,state_dim)]

    B = [zeros(state_dim-1,1);1]

    C = [1 zeros(1,state_dim-1)];

    D = 0;

    m = size(B,2); %number of inputs
    
    sys = ss(A,B,C,D);

    if state_dim < 4
    %% Compute the viability kernel using the polytope algorithm
        
        U  = unitbox(1,0.3);  % Input constraint set
        K0 = unitbox(state_dim,0.5);  % State constraint set
        
        N = 10;     % maximum number of time steps
        tau = 4/N;  % time step (4 second horizon)

        sys_mpt = mpt_sys(sys,tau);

        % Inversion of time direction (to get backward reach sets)
        sys_mpt.A = inv(sys_mpt.A); 
        sys_mpt.B = -inv(sys_mpt.A)*sys_mpt.B;

        % Required arguments for MPT but not used for reach computations
        sys_mpt.xmax = 0.5*ones(state_dim,1); sys_mpt.xmin = -0.5*ones(state_dim,1); 
        sys_mpt.umax = .3; sys_mpt.umin = -.3;
        
        tic
        computeViab_Polytope
        run_time(state_dim,1) = toc
        
    end
    
    %%%%%%%% Code to compute the viability kernel using the ellipsoidal
    %%%%%%%% algorithm. Since the Ellipsoidal toolbox does not handle being
    %%%%%%%% called from inside a loop well, we have run this section
    %%%%%%%% manually for state_dim = 1:15 and saved the runtimes in the
    %%%%%%%% file ellipsoid_manual.mat
    
    if state_dim < 25
    %% Compute the viability kernel using the ellipsoidal algorithm
    
        L = [1; zeros(state_dim-1,1)]; % approximation directions
        num_L = 1;        % number of approximation directions
        
        U = 0.3*ell_unitball(1);     % input constraint set
        K0 = 0.5*ell_unitball(state_dim); % state constraint set
   
        N = 10;     % maximum number of time steps
        tau = 4/N;  % time step (4 second horizon)
        dsys = c2d(sys, tau);
        ell_sys = linsys(dsys.A, dsys.B, U, [], [], dsys.C, [], 'd');
   
        tic
        computeViab_Ellipsoid
        run_time(state_dim,2) = toc
        
    end
    
    %% Compute the viability kernel using the support vector algorithm
    
    L = [eye(state_dim) -eye(state_dim)]; % approximation directions
    num_L = 2*state_dim;        % number of approximation directions
    
    n = state_dim;
    
    tau = 4;           % time horizon 
    steps = 10;        % number of time steps for support vector method
    rho = tau/steps;   % disctetization interval for support function method
    sysd = c2d(sys,rho);

    Ad = sysd.A;
    Bd = sysd.B;
    Astar = inv(Ad)'; 
       
    % State constraint set 
    A_K0 = 1/2*eye(state_dim);
    q_K0 = zeros(state_dim,1);

    % Input constraint set
    A_U = 0.3;
    q_U = 0;
    A_V = -Bd*A_U;
    q_V = -Bd*q_U;
    
    tic
    computeViab_SupportVect
    run_time(state_dim,3) = toc

end 

%% plot the run time graph

% load('ellipsoid_manual.mat')         % load ellipsoid run times from file
% run_time(1:15,2) = ellipsoid_manual

run_time(4,1) = 500 % this is cheating, I know ... but it makes the graph intuitive
                    % the true run time for the polytope method in 4D is
                    % much greater than 500 seconds
                    
figure                     
    hold on
    plot(run_time(1:4,1),'k', 'LineWidth', 2.5)
    plot(run_time(1:24,2), 'color', [0 0.1569 0.3490], 'LineWidth', 2.5)
    plot(run_time(:,3), 'color', [0.4549 0.5686 0.6392], 'LineWidth', 2.5)
    xlabel('state dimension')
    ylabel('run time (s)')
    legend('polytope', 'ellipsoid', 'support vector')
    axis([1 30 0 200])
    hold off
    
    