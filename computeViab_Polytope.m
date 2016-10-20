%% Computes exact viability kernel using polytope method
%   viability kernel for time horizon N is given as polytope KN
%
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires Multi-Parametric Toolbox

Knminus1 = polytope([]);
Kn = K0;
n = 0;

while n <= N
    if  isempty(Kn)
        KN = polytope([]);
        disp('Empty K0')
        break
    end
    if  eq(Kn,Knminus1)
        KN = Kn;
        disp(['Viability kernel converged after ' num2str(n) ' time steps' ])
        %disp(num2str(n))
        break
    end
    L = mpt_reachSets(sys_mpt, Kn, U, 1); % The one-step backward reach set
    Knminus1 = Kn;
    Kn = L & K0;
    n = n+1;
end
KN = Kn;
