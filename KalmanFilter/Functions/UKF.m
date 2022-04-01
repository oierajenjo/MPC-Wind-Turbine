function [xk,P,e] = UKF(f,h,Q,R,xk,y_me,u_b,Lk,Yk,N,P0,Ts,v,n)
%% Kalman variables
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(Lk+kappa) - Lk;
n_sigma_p = 2*Lk + 1; % Number of sigma points
wm = ones(n_sigma_p,1)*1/(2*(Lk+lambda)); % Weight for transformed mean
wc = wm; % Weight for transformed covariance
wm(1) = lambda/(lambda+Lk);
wc(1) = lambda/(lambda+Lk) + 1 - alpha^2 + beta;

%% Unscented Kalman Filter
P = zeros(N,Lk,Lk);
P = P0; % Set initial error covariance
e = zeros(Yk,N);
for k = 1:N-1
    % Step 1: Generate the sigma-points
    try 
        sP = chol(P,'lower');% Calculate square root of error covariance
%       norm(sP*sP' - P)
%         disp('Matrix is symmetric positive definite.');
    catch ME
        disp('Matrix is not symmetric positive definite');
        k
%         % 1) Compute the LDL decomposition: A lower-triangular matrix L, 
%         %    a block-diagonal matrix D (1-by-1 and 2-by-2 blocks),
%         %    and a permutation matrix P, such that A is P*L*D*L'*P'
%         [l, d, p] = ldl(P);
%         % 2) Compute the eigenvalue decomposition, set its negative eigenvalues to zero,
%         % and use QR to get two triangular factors for this modified matrix:
%         [U, d] = eig(P, 'vector'); % A == U*diag(d)*U'
%         d(d < 0) = 0; % Set negative eigenvalues of A to zero
%         [~, ~] = qr(diag(sqrt(d))*U'); % Q*R == sqrt(D)*U', so A == U*diag(d)*U'
%                                     %                         == R'*Q'*Q*R == R'*R
%         % In this case, check that d(d<0) are all nearly zero.
%         % 3) Ugly but quick: Just add a small number to A's diagonal before calling Cholesky
%         sP = chol(P + 1e-7*eye(size(P)),'lower');
        break
    end
     
    % chi_p = "chi previous" = chi(k-1) % Untransformed sigma points
    chi_p = [xk(:,k), xk(:,k)*ones(1,Lk)+sqrt(Lk+lambda)*sP, ...
        xk(:,k)*ones(1,Lk)-sqrt(Lk+lambda)*sP]; % Untransformed sigma points
    
    % Step 2: Prediction Transformation
    % Propagate each sigma-point through prediction
    % chi_m = "chi minus" = chi(k|k-1)
    chi_m = zeros(Lk,n_sigma_p); % Transformed sigma points
    for j=1:n_sigma_p
        % Runge-Kutta 4th order method
        k_1 = f(chi_p(:,j),u_b(:,k));
        k_2 = f(chi_p(:,j)+0.5*Ts*k_1,u_b(:,k)+0.5*Ts);
        k_3 = f(chi_p(:,j)+0.5*Ts*k_2,u_b(:,k)+0.5*Ts);
        k_4 = f(chi_p(:,j)+Ts*k_3,u_b(:,k)+Ts);
        chi_m(:,j) = chi_p(:,j) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n(chi_p(:,j)); 
%         chi_m(:,j) = chi_p(:,j) + Ts*f(chi_p(:,j),u_b(:,k));
    end
    
    x_m = chi_m*wm; % Calculate mean of predicted state
    % Calculate covariance of predicted state
    P_m = Q(xk(:,k)); % A priori covariance estimate
    for j = 1:n_sigma_p
        P_m = P_m + wc(j)*(chi_m(:,j) - x_m)*(chi_m(:,j) - x_m)';
    end
    
    % Step 3: Observation Transformation
    % Propagate each sigma-point through observation
    % Initial velocity will be considered as 0, as we need it for
    % obtaining the acceleration
    psi_m = zeros(Yk,n_sigma_p);
    for j=1:n_sigma_p
        psi_m(:,j) = h(chi_m(:,j)) + v(:,k);
    end
    y_m = psi_m*wm; % Calculate mean of predicted output
    
    % Calculate covariance of predicted output
    % and cross-covariance between state and output
    Pyy = R;
    %Pxy = zeros(L,2);
    Pxy = 0;
    for j = 1:n_sigma_p
        Pyy = Pyy + wc(j)*(psi_m(:,j) - y_m)*(psi_m(:,j) - y_m)';
        Pxy = Pxy + wc(j)*(chi_m(:,j) - x_m)*(psi_m(:,j) - y_m)';
    end
    
    % Step 4: Measurement Update
    K = Pxy/Pyy; % Calculate Kalman gain
    e(:,k+1) = y_me(:,k+1) - y_m;
    xk(:,k+1) = x_m + K*e(:,k+1); % Update state estimate
    P = P_m - K*Pyy*K'; % Update covariance estimate
end
end

