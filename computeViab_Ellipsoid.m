%% Computes exact viability kernel using ellipsoidal method
%   various under-approximations of the viability kernel are given in the
%   array of ellipsoids KN
%
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires Ellipsoidal Toolbox

Knminus1 = polytope([]);
Kn = K0;
n = 0;

for j = 1:num_L
    fprintf(['\n =======Computing ellipsoid in direction ' num2str(j) ' of ' num2str(num_L) '=======\n'])
    ell = L(:,j);
    n = 0;
    Knminus1 = ellipsoid;
    Kn = K0;
    while n <= N
        if  isempty(Kn)
            KN(j) = ellipsoid;
            disp('Empty K_N')
            break
        end
        if  eq(Kn,Knminus1)
            KN(j) = Kn;
            disp(['Viability kernel converged after ' num2str(n) ' time steps' ])
            break
        end
        rs = reach(ell_sys, Kn, ell, [1 0], options);  % compute reach set
        L_int = get_ia(rs);
        Knminus1 = Kn;
        Kn = ellintersection_ia([L_int  K0]);          % intersect with state constraints
        n = n+1;
    end
    KN(j) = Kn;
end
