%% Approximation of the viability kernel using support functions and support vectors
%     compared with level set approximation
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires CVX, Ellipsoidal Toolbox, Multi-Parametric
%     Toolbox


clear all; clc
cvx_quiet(true);

%% Continuous-time double integrator

A = [0 1;
     0 0];
 
B =  [ 0 ;
       1  ]; 
   
C = [1 0];

D = 0;
   
sys = ss(A,B,C,D);
n = size(sys.A,2);
m = size(sys.B,2);


%% Discretization

tau = 4;           % time horizon 
steps = 20;        % number of time steps for support vector method
rho = tau/steps;   % disctetization interval for support function method

sysd = c2d(sys,rho);

Ad = sysd.A;
Bd = sysd.B;
Astar = inv(Ad)'; 
        

%% State constraint set 
%    image of unit ball under x |-> (A_K0 x + q_K0)

A_K0 = 0.5*eye(n);
q_K0 = zeros(n,1);


%% Input constraint set
%    image of unit ball under x |-> (A_V x + q_V) in v-space

A_U = 0.3;
q_U = 0;

A_V = -Bd*A_U;
q_V = -Bd*q_U;


%% Setup approximation directions 

num_L = 20; % number of support vectors

theta = linspace(0,2*pi, num_L+1);
for j=1:num_L
L(:,j) = [cos(theta(j)) -sin(theta(j));
          sin(theta(j))  cos(theta(j))] * [1;0];
end


%% Run support vector algorithm 

tic

computeViab_SupportVect

time_6b = toc


%% compute discrete viability kernel via gridding
% commented code does this, but takes a long time to run
% instead of running it each time we simply load the 
% results from 'grid_viab.mat'

[X,Y] = meshgrid(-0.5:0.001:0.5,-0.5:.001:0.5);

% Z = -ones(size(X));

% u = 0.3;

%for i=1:size(X,1)
 %   for j=1:size(X,2)
  %      x = [X(i,j);
   %          Y(i,j)];
    %    count = 0;
     %   if (Y(i,j) < 0)
      %      while(Z(i,j) == -1)
       %           count = count+1;
        %        if( norm(x) >= 0.5 )
         %           Z(i,j) = 0;
          %      end
           %     if( x(2) >= 0 || count >= steps)
            %        Z(i,j) = 1;
             %   end
              %  x = Ad * x - Bd * sign(x(2)) * u;
%            end
 %       end
  %      if (Y(i,j) >= 0)
   %         while(Z(i,j) == -1)
    %             count = count+1;
     %           if( norm(x) >= 0.5 )
      %              Z(i,j) = 0;
       %         end
        %        if( x(2) <= 0 || count >= steps)
         %           Z(i,j) = 1;
          %      end
           %     x = Ad * x - Bd * sign(x(2)) * u;
            %end
%        end
%    end
%end

load('grid_viab.mat')


%% Plot the results

Options1.color=[0 0.1569 0.3490]; % UBC blue
Options2.color=[0.4549 0.5686 0.6392]; % UBC grey

P_over = polytope(A_viab,b_viab);
P_under = polytope(v);

figure
plot(P_over, Options2)
hold on 
contour(X,Y,Z,[1 1],'LineWidth',2,'fill','on')
colormap([0 0 0])
plot(P_under, Options1)
axis([-.6,.6,-.6,.6])
set(gca,'DataAspectRatio',[1 1 1])  

    
    