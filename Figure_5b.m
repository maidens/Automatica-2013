%% Approximation of the viability kernel using ellipsoids 
%    (compared with level set appooximation)
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires Ellipsoidal Toolbox 

clear all; clc


%% Continuous-time double integrator

A = [0 1;
     0 0];
 
B =  [ 0 ;
       1  ]; 
   
C = [1 0];

D = 0;
   
sys = ss(A,B,C,D);

U = 0.3*ell_unitball(1);     % input constraint set

K0 = 0.5*ell_unitball(2); % state constraint set


%% Discretization

N = 40;     % maximum number of time steps
tau = 4/N;  % time step (4 second horizon)

dsys = c2d(sys, tau);

ell_sys = linsys(dsys.A, dsys.B, U, [], [], dsys.C, [], 'd');


%% Ellipsoidal Toolbox options

global ellOptions;
ellOptions.verbose = 0;
options.approximation = 1;  %0: external, 1:internal, 2:both


%% Define Ellipsoidal Toolbox approximation directions

num_L = 20;    % number of approximation directions

theta = linspace(0, pi, num_L+1);
for j=1:num_L
L(:,j) = [cos(theta(j)) -sin(theta(j));
          sin(theta(j))  cos(theta(j))] * [1;0];
end


%% Run viability kernel approximation algorithm

tic

computeViab_Ellipsoid

time_5b = toc


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
        %        if( norm(x) >= N )
         %           Z(i,j) = 0;
          %      end
           %     if( x(2) >= 0 || count >= N)
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
        %        if( x(2) <= 0 || count >= N)
         %           Z(i,j) = 1;
          %      end
           %     x = Ad * x - Bd * sign(x(2)) * u;
            %end
%        end
%    end
%end

load('grid_viab.mat')

    
%% Plot the results

Options1.color=[0 0.1569 0.3490];      % UBC blue
Options1.fill = 1;
Options2.color=[0.4549 0.5686 0.6392]; % UBC grey
Options2.fill = 1;

figure
    plot(K0,Options2)           % plot state constraint set
    hold on
    contour(X,Y,Z,[1 1],'LineWidth',2,'fill','on')
    colormap([0 0 0])
    for l=1:num_L
        plot(KN(l),Options1)    % plot ellipsoidal underapproximations of viability kernel
    end
    axis([-.6,.6,-.6,.6])
    set(gca,'DataAspectRatio',[1 1 1])
    
   



