%% Computes approximation of viability kernel using support vector method
%   for heat equation example 
%
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires computeViab_SupportVect, CVX, Ellipsoidal Toolbox,
%    Multi-Parametric Toolbox

clear all; clc
cvx_quiet(true)


%% continuous-time state matrices

n = 20
alpha = 8;
dx = 1;

A = -2.014*eye(n);
A(1,1) = -1;
A(n,n) = -1;
for i=1:n-1
    A(i,i+1) = 1;
    A(i+1,i) = 1;
end
A = alpha/dx*A;

B = zeros(n,1);
B(1) = 1;

C = zeros(1,n);

D = 0;

sys = ss(A,B,C,D);

n = size(sys.A,2);
m = size(sys.B,2);

%% discretization

tau = 5;          % time horizon 
steps = 20;        % number of time steps for support function method
rho = tau/steps;   % disctetization interval for support function method

sysd = c2d(sys,rho);

Ad = sysd.A;
Bd = sysd.B;
Astar = inv(Ad)';


%% state constraint set is image of unit ball under x |-> (A_K0 x + q_K0)
%    units of mcg/mL

A_K0 = diag(10*ones(n,1)); 
q_K0 = 10*ones(n,1);


%% input constraint set is [0, u_max] (in u-space)
%    (i.e. image of unit ball under x |-> (A_V x + q_V) in v-space)

A_U = 5*diag(ones(m,1));

A_V = -Bd*A_U;
q_V = -Bd*A_U;


%% set directions for support vector method

% box over-approximation / dual box under-approximation
% L = [eye(n) -eye(n)];

% evenly
num_L = 12; % number of supporting hyperplanes
theta = linspace(0,2*pi, num_L+1);
theta = theta(1:end-1);
for j=1:num_L
L1(:,j) = [cos(theta(j)) -sin(theta(j));
          sin(theta(j))  cos(theta(j))] * [1;0];
end

L = zeros(n,n/2*num_L);
for k=1:n/2
    L(2*k-1:2*k,(k-1)*num_L+1:k*num_L)=L1;
end


L = A_K0\L;
num_L = size(L,2);

%% Run Algorithm 4 (support vector method) to compute overapproximating polytope
%   (A_viab, b_viab) and underapproxiamting polytope v for viability kernel

tic

computeViab_SupportVect

time_9 = toc;


%% Construct geometric objects using Ellipsoidal and Multi-Parametric 
%    Toolboxes

P = polytope(A_viab,b_viab);

for i=1:(n/2)
    U(i) = A_K0(2*i-1:2*i,2*i-1:2*i)*ell_unitball(2)+q_K0(2*i-1:2*i);
    v_proj = v(:,2*i-1:2*i);
    P_under(i) = polytope(v_proj); 
    P_over(i) = projection(P,2*i-1:2*i);
end

U_star = U(1);

v_proj1 = [v(:,1) v(:,20)];
P_under_star1 = polytope(v_proj1);
P_over_star1 = projection(P,[1 20]);

v_proj2 = [v(:,1) v(:,10)];
P_under_star2 = polytope(v_proj2);
P_over_star2 = projection(P,[1 10]);

%% Plot the above objects

Options1.color=[0 0.1569 0.3490];        % UBC blue
Options2.color=[0.4549 0.5686 0.6392];   % UBC grey

figure 
   subplot(2,2,1)
   hold on
        plot(P_over(1), Options2)
        plot(P_under(1),Options1)
        plot(U(1))
        xlabel(strcat('\xi_{1}'))
        ylabel(strcat('\xi_{2}'))
    hold off
    
   subplot(2,2,2)
   hold on
        plot(P_over(10), Options2)
        plot(P_under(10),Options1)
        plot(U(10))
        xlabel(strcat('\xi_{19}'))
        ylabel(strcat('\xi_{20}'))
    hold off

    
   subplot(2,2,3)
   hold on
        plot(P_over_star2, Options2)
        plot(P_under_star2,Options1)
        plot(U_star)
        xlabel(strcat('\xi_{1}'))
        ylabel(strcat('\xi_{10}'))
    hold off

    
   subplot(2,2,4)
   hold on
        plot(P_over_star1, Options2)
        plot(P_under_star1,Options1)
        plot(U_star)
        xlabel(strcat('\xi_{1}'))
        ylabel(strcat('\xi_{20}'))
    hold off


 

