%% Computes approximation of viability kernel using support vector method
%   under-approximation is given as the V-polytope Conv(v)
%   over-approximation is given as the H-polytope {x : A_viab*x < b_viab}
%
%  John Maidens 2012 (jmaidens@ece.ubc.ca)
%  Requires CVX

A_viab = zeros(num_L,n);
b_viab = zeros(num_L,1);
v = zeros(num_L,n);

for j=1:num_L  % loop over each direction in the set L 
    fprintf(['\n =======Computing support vector in direction ' num2str(j) ' of ' num2str(num_L) '=======\n'])
    l = L(:,j);
    clearvars xi
    
    
    %% Minimize the function xi defined in (21) using CVX
    %     note: support functions rho_K0 and rho_V given for ellipsoidal
    %     input and state constraints; 
    %     easily adapted to other support functions 
    
    cvx_begin sdp
    cvx_precision best
    cvx_solver sdpt3
        variable w(n*steps)
        for k=1:steps-1
            xi(k) = norm(A_K0'*(Astar*w((k-1)*n+1:k*n) - w(k*n+1:(k+1)*n)),2)   +   q_K0'*(Astar*w((k-1)*n+1:k*n) - w(k*n+1:(k+1)*n))   +   norm(A_V'*(Astar*w((k-1)*n+1:k*n)),2)  +   q_V'*(Astar*w((k-1)*n+1:k*n));
        end
       minimize(sum(xi) + norm(A_K0'*(l-w(1:n)),2) + q_K0'*(l-w(1:n)) + norm(A_V'*(Astar*w(end-n+1:end)),2) + q_V'*(Astar*w(end-n+1:end)) + norm(A_K0'*(Astar*w(end-n+1:end)),2) + q_K0'*(Astar*w(end-n+1:end)))
    cvx_end

    
    %% Store optimization value for over-approximation of the viability kernel
    
    A_viab(j,:) = l';        % store direction l
    b_viab(j) = cvx_optval;  % store value of the support function rho(l)

    
    %% Use the value of the optimization parameter w to compute and store a support vector v(l) 
    %     note: this implementation works only for ellipsoidal 
    %     input and state constraints
    
    if norm(l-w(1:n)) > 2* cvx_slvtol
        % store support vector v(l) on the boundary of the viability constraint set
        v(j,:) = q_K0' + (l-w(1:n))'*(A_K0*A_K0')/norm(A_K0'*(l-w(1:n)),2);
    else      
        k = 1;
        iter = (q_V' + (Astar*w(1:n))'*(A_V*A_V')/norm(A_V'*(Astar*w(1:n)),2))*Astar;
        while norm(Astar*w((k-1)*n+1:k*n) - w((k)*n+1:(k+1)*n),2) < 2*cvx_slvtol  && (k < tau) 
           k = k+1;
           iter = iter +  (q_V' + (Astar*w((k-1)*n+1:k*n))'*(A_V*A_V')/norm(A_V'*(Astar*w((k-1)*n+1:k*n)),2)) * Astar^k;
        end
        % store support vector v(l) not on the boundary of the viability constraint set
        v(j,:) = iter + (  q_K0' + (Astar*w((k-1)*n+1:k*n) - w((k)*n+1:(k+1)*n))'*(A_K0*A_K0')/norm(A_K0'*(Astar*w((k-1)*n+1:k*n)-w((k)*n+1:(k+1)*n)),2)  )*Astar^k;
        if k == tau
            warning(strcat('No value set for direction l', num2str(j)))
        end        
    end
   
    
end
