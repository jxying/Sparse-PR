function out =  solver_newton_projection(y_abs,A,s,opts)
% Sparse phase retrieval using sparse spectral initialization method in the first stage
% and the Newton projection algorithm in the second stage.

% Inputs:
%    'y_abs'    : phaseless measurements
%    'A'        : sensing matrix, where each column is an sensing vector;
%    'opts'     : stopping criterion options, passed as a structure possible field name values:             
%                 'maxiter'   : maximum number of iterations;  
%                 'tol'       : tolerance, the algorithm stops if dist(x_k+1 - x_k) / ||x_k|| < tol 
%                               if the fround true signal is not available for evaluating termination condition;  
%                 'tolpcg'    : tolerence of stopping criterion for the preconditioned congurate gradient algorithm;
%                 'true'      : the ground truth signal, input for evaluating termination condition; 
%                 'trunction' : "1" if we use truncated measurements, 
%                               "0" otherwise (aligns with the setting in the paper)
%                 'svd'       : 'svd'for a more accurate solution in svd during the first stage but slower in high dimensions, 
%                               'power' for approximate solution but faster;
%                 'eta'       : stepsize in hard thresholding, "0.95" by default;
%                 'tau'       : stepsize in Newton projection, "1" by default (aligns with the setting in the paper for Quacdratic convergence)
%                 'display'   : "1" if display the iteration results, "0" otherwise.
%
% 
% Outputs:
%     'out' : outputs as a structure
%             possible field name values:
%                'x'             : estimated signal; 
%                'time'          : run time cost before algorithm stopped
%                'timeiter'      : store the cpu time for each iteration
%                'iterate'       : the number of iterations cost before algorithm stopped 
%                'relerror'      : relative error  of the estimated signal, i.e., dist(x_k - x^*) / \| x^* \|,
%                                    exists if the ground truth signal is available;   
%                'relerror_iter' : relative error stored for each iteration, 
%                                    exists if the ground truth signal is available;
%                'converge'      : "1": algorithm converges; "0" algorithm does not converge.


%% Set default parameters
if nargin < 4
   opts.maxiter  = 1e2;
   opts.tol      = 1e-5;
   opts.tolpcg   = 1e-4;
end

if isfield(opts,'maxiter')
    maxiter = opts.maxiter;
else
    maxiter = 100;
end

if isfield(opts,'tol')
    tol = opts.tol; 
else
    tol = 1e-5;
end

if isfield(opts,'tol_pcg')
    tolpcg = opts.tolpcg; 
else
    tolpcg = 1e-4;
end

if isfield(opts,'true')
    z = opts.true; 
    flag_true = 1;
    if isfield(opts,'toltrue')
       toltrue = opts.toltrue;
    else
       toltrue = 1e-3;
    end
else
    flag_true = 0;
end

if isfield(opts,'trunction')
    flag_trunction = opts.trunction;
else
    flag_trunction = 0;
end

if isfield(opts,'svd')
    option_svd = opts.svd;
else
    option_svd ='power';
end

if isfield(opts,'eta')
    eta = opts.eta;
else
    eta = 0.95;
end

if isfield(opts,'tau')
    tau = opts.tau; 
else
    tau = 1;
end

if isfield(opts,'display')
    display     = opts.display;
else
    display     = 1;
end

%% Initialize parameters
converge  = 0;
succ      = 0;
[m,n]     = size(A);
xo        = zeros(n,1);
erroriter = zeros(maxiter,1); % relative error stored for each iteration
timeiter  = zeros(maxiter,1); % cpu time stored for each iteration
MShat     = zeros(s);         
y_abs2    = y_abs.^2;         % measurements
phi_sq    = sum(y_abs2)/m;
phi       = sqrt(phi_sq);     

%% First Stage
Marg     = ((y_abs2)'*(A.^2))'/m; 
[Mg,MgS] = sort(Marg,'descend');
Shat     = sort(MgS(1:s));        % estimate support of the signal
AShat    = A(:,Shat);             

if flag_trunction                 % using truncated measurements
   card_Marg = ceil(m/6);
   for i=1:m
       M_eval(i) = y_abs(i)/norm(AShat(i,:));
   end 
   [~,MmS] = sort(M_eval,'descend');
   Io = MmS(1:card_Marg);   
else                              % using untruncated measurements 
   card_Marg = m;
   Io = 1:card_Marg;
end

for i = 1:card_Marg
    ii = Io(i);
    MShat = MShat + (y_abs2(ii))*AShat(ii,:)'*AShat(ii,:); 
end

switch option_svd
    case 'svd'
        [u,sigma,v] = svd(MShat);
        v1 = u(:,1);     %top singular vector
    case 'power'
        v1 = svd_power(MShat);
end

v = zeros(n,1);
v(Shat,1) = v1;
x_init = phi*v;         % rescaling the initial estimate
x = x_init;

%% Second Stage
t0 = tic;               % initialize CPU time for the second stage
Tx      = find(x); 
Ax      = A(:,Tx)*x(Tx);
Ax3     = Ax.^3;
Axb3    = Ax3-y_abs2.*Ax;
g       = A'*Axb3/(m*m);
T0      = Tx;

for t=1:maxiter 
    
    x0      = x;
    p       = sign(Ax);
    y_sign  = y_abs.*p;
    grad_i  = A'*(Ax-y_sign)/m;

    xtg     = x0-eta*grad_i;
    [~,T]   = maxk(abs(xtg),s); 
    
    TTc     = setdiff(T0,T);
    T0      = T;
    gT      = g(T);
    
    Ax3b    = 3*Ax.^2-y_abs2;
    AT      = A(:,T); 
    H       = AT'*(Ax3b.*AT)/(m*m); 

    D       = AT'*(Ax3b.*A(:,TTc))/(m*m);
    [dT,~]  =  pcg(H, D*x0(TTc)-gT, tolpcg);
    
    x       = xo;
    x(T)    = x0(T) + tau*dT;

    Tx      = find(x); 
    Ax      = A(:,Tx)*x(Tx);
    Ax3     = Ax.^3;
    Axb3    = Ax3-y_abs2.*Ax;
    g       = A'*Axb3/(m*m);

    timepoint    = toc(t0);
    timeiter(t)  = timepoint;

    if flag_true
       err = compute_error(x,z);
       erroriter(t) = err;
       if display
             disp(['iter: ', num2str(t), '   ', 'relative error: ', num2str(err,11), '   ', 'cpu time: ', num2str(timepoint) ]);
       end
       if err < toltrue
           converge = 1;
           break;
       end
       if err > 0.5
           succ = succ + 1;
       else
           succ = 0;
       end

       if succ >= 10
           disp('signal recovery failed: please increase measurements or regenerate sensing matrix');
           break
       end

    else
       err = compute_error(x,x0);
       if display
             disp(['iter: ', num2str(t), '   ', 'difference between successive iterations: ', num2str(err,11), '   ', 'cpu time: ', num2str(timepoint) ]);
       end
       if err < tol
           converge = 1;
           break;
       end
    end

end

out.x                = x;
out.time             = timepoint;
out.timeiter         = timeiter(1:t);
out.iterate          = t;
out.converge         = converge;
if flag_true
   out.relerror      = err;
   out.relerror_iter = erroriter(1:t);
end

end











