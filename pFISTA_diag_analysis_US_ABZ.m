function [ x ] = pFISTA_diag_analysis_US_ABZ( y, H, Params )
% PFISTA_DIAG_ANALYSIS_US - Performs the Fast Proximal Gradient method on the problem
% 
%              min_{x} \lambda|| W^*x ||_1 + 0.5|| Y - AXA^H ||_F^2
%   
% with X being diagonal with a sparse diagonal x. A is assumed to have the following 
% specific structure: A = H*kron(Fp, Fp), with H being a diagonal matrix  and Fp 
% is a partial Fourier matrix. This structure leads to a very efficient implementation 
% (in terms of run-time and memory usage) is using fft and ifft operations. 
%
% It is also assumed that x is non-negative and real (which is true for a correlation matrix, for instance)
%
% Syntax:
% -------
% [ x ] = pFISTA_diag_analysis_US( Y, H, Params )
%
% Inputs:
% -------
% Y      - Input M^2 X M^2 observations matrix (Generaly complex)
% H      - M^2 X M^2 diagonal matrix (should be in sparse format)
% Params - Additional algorithmic parameters
%          Beta      : Regularization parameter (0,1) - best to take very close to 1
%          L0        : Lipschitz constant (norm(kron(A, A), 2))
%          LambdaBar : Regularization parameter. Should be very small, e.g. 1e-7
%          Lambda    : Sparsity regularization parameter. Value depends on the problem
%          N         : Reconstructed Xr is of size N^2 X N^2
%          isPos     : 1 if X entries are supposed to be non negative
%          IterMax   : Maximum number of iterations
%          LargeScale: 1 if N > 500^2 (roughly). Computation slows dramatically, but can handle very large datasets
%
% Output:
% -------
% x      - An estimate of the sparse N X 1 vector.
%
% Written by Oren Solomon, Technion I.I.T. Ver 1. 07-07-2016
%

global VERBOSE

%% Initializations
% -----------------------------------------------------------------------------------------------------------------------
% Paramters
Calc_S_flag = 0;
beta       = Params.Beta;
L          = Params.L0;                % Lipschitz constant
lambda_bar = Params.LambdaBar;
lambda     = Params.Lambda;            % l1 regularization parameters
N          = Params.N;                 % Length of x
mu         = Params.mu;

% Init
F_stack    = zeros(Params.IterMax, 1); % Accumulate function values

t          = 1;

% Memory allocation
x          = zeros(N, 1);
z          = x;

% Determine type of analysis operator
switch lower(Params.AnalysisType)
    case 'wave1d' % Wavelet transform - based on the Rice Wavelet Toolbox, Version 3.0, Dec 2002 - 1D
        WaveD     = 2;                                      % Wavelet depth hierarchy
        WaveScale = daubcqf(8, 'mid');                     % Type of wavelet filter - Daubechies wavelets
        D         = @(x) midwt(x, WaveScale, WaveD);        % Synthesis operator
        Dt        = @(x) mdwt(x, WaveScale, WaveD);         % Analysis operator
        
    case 'wave2d' % Wavelet transform - based on the Rice Wavelet Toolbox, Version 3.0, Dec 2002 - 2D
        WaveD     = 1;%1;                                      % Wavelet depth hierarchy - Daubechies wavelets
        WaveScale = daubcqf(16, 'mid'); %64                    % Type of wavelet filter
        D         = @(x) vec(midwt(reshape(x, sqrt(N), sqrt(N)), WaveScale, WaveD));  % Synthesis operator
        Dt        = @(x) vec(mdwt(reshape(x, sqrt(N), sqrt(N)), WaveScale, WaveD));   % Analysis operator
        
    case 'fwht' % Fast Welch Hadamard transform
        D  = @fwht;                                         % Synthesis operator
        Dt = @ifwht;                                        % Analysis operator
        
    case 'dct1d'  % DCT - 1D
        D  = @dct;                                          % Synthesis operator
        Dt = @idct;                                         % Analysis operator
        
        % Dnorm = norm(dctmtx(N), 2);
    case 'dct2d'  % DCT - 2D
        D  = @(x) vec(dct2(reshape(x, sqrt(N), sqrt(N))));  % Synthesis operator
        Dt = @(x) vec(idct2(reshape(x, sqrt(N), sqrt(N)))); % Analysis operator
        
        % Dnorm = norm(dctmtx(N), 2);
    otherwise
        error('pFISTA_diag_analysis: Unknown analysis operator.');
end

Dnorm = 0.1; 0.01;

if VERBOSE; disp('pFISTA_diag: Calculating the Lipschitz constant...'); end;

% Calculate S
if VERBOSE; fprintf('S_calc: ');t2 = tic; end;
if Calc_S_flag
    S = S_calc(H, N); % ???
else
    switch Params.N 
        case (8^2)*(64^2)
            load('Smat_64_to_512.mat');
        case (16^2)*(64^2)
            load('Smat_64_to_1024.mat');
        otherwise
            load('Smat_64_to_512.mat');
    end
end
if VERBOSE; toc(t2); end;

if VERBOSE; fprintf('v_calc: ');t3 = tic; end;
v = v_calc(y, H, N); 
if VERBOSE; toc(t3); end;

% Lipschitz constant
if isempty(L)
    L = max(max(S)) + (Dnorm^2)/mu;
end

lambda_mu = lambda*mu;

if VERBOSE; disp('done.'); end;
Titers = 0; 

%% Iterations
% -----------------------------------------------------------------------------------------------------------------------
if VERBOSE; disp('pFISTA_diag: Running iterations...'); end;
for kk = 1:Params.IterMax
    tic;
    if VERBOSE && mod(kk,50)==0; fprintf(['Iteration #' num2str(kk) ': ']); end;
      
    % Gradient step: G = Z - (1/L)*( |A^H*A|^2*x - v + Grad_g )
    g = z - (1/L)*(Grad_f(z, S, v, N) + Grad_g(x, D, Dt, lambda_mu, mu));

    % Apply monotonicity (Monotone FISTA)
    x_prev = x;                                                        % x_prev = x_{k-1}
    
    if Params.mon
        Fval_tmp = FuncValue(x, Y, Dt, mu, Lambda);
        if Fval_tmp > FuncValue(g, Y, Dt, mu, Lambda)     % x = argmin{H_mu(x): z_k, x_{k-1}}
            x = g;
        end
        
        % Update function value stack
        F_stack(kk) = Fval_tmp;
    else
        x = g;
    end
    
    % Projection onto the non-negative orthant, only for a non-negative constraint
    if Params.NonNegOrth == 1
        x(x < 0) = 0;
    end
    
    % Parameter updates for next iteration
    t_prev = t;
    t = 0.5*(1 + sqrt(4*t^2 + 1));
    
    % Z update
    z = x + ((t_prev - 1)/t)*(x - x_prev) + (t_prev/t)*(g - x);
    
    lambda = max(beta*lambda, lambda_bar);
    
    Titers = Titers + toc; 
    if VERBOSE && mod(kk,50)==0; 
        disp(['Elapsed time is ' num2str(Titers) ' seconds.'])
        Titers = 0; 
    end;
end
if VERBOSE; disp('Done pFISTA_diag.'); end;

%% Auxiliary functions
% -----------------------------------------------------------------------------------------------------------------------
%% Soft thresholding
function y = Soft(z, alpha)
y = sign(z).*max(abs(z) - alpha, 0);

%% Compute fanction value
function f = FuncValue(x, Y, Dt, mu, lambda)
AXA_H = LAH(diag(x), H);                        % AXA_H = A*X
AXA_H = LAH(AXA_H', H);                         % AXA_H = A*AXA_H'
f_val = 0.5*norm(Y - AXA_H', 'fro')^2;
g_val = lambda*sum(Huber(Dt(x), lambda*mu));

f = f_val + g_val;

% %% Implementation of a_i^H*p
% function a = lkfft(ii, Pv, Nsqrt, Msqrt)
% % Step 0: Determine indices
% ki = floor((ii - 1)/Nsqrt) + 1;
% li = mod(ii, Nsqrt) + Nsqrt*(mod(ii, Nsqrt) == 0);
% 
% % Step 1: Convert to M X M matrix
% t = reshape(Pv, Msqrt, Msqrt);
% 
% % Step 2: Q is an N X M matrix
% Q = Nsqrt*ifft(t, Nsqrt); % N*ifft(t, Nsqrt); Why Nsqrt^2 and not Nsqrt ?
% 
% % Step 3: Output is an N X 1 vector
% q2 = fft(Q(li, :)', Nsqrt);
% 
% % Output
% a = q2(ki);

function v = v_calc(y, H, N)
v = LAH_H( y, H, N );

%% Calculate A1 which is needed for the gradient calculation 
function S = S_calc(H, N)
% Step 0
H2  = diag( abs(diag(H)).^2 );

% Step 1: M^2 X 1 vector q
q  = ctranspose( LAH_I(H2, N) );

% Step 2: N^2 X 1 vector A1
A1 = LAH_H(q, 1, N);

% Step 3: Calculate eigenvalues sqrt(N) X sqrt(N) matrix (N eigenvalues)
S = fft2(reshape(A1, sqrt(N), sqrt(N)));        % 1/N ? - does not seem to affect reconstruction performance

%% Calculate the gradient step: \nabla f(x) = |A^H*A|^2*x - v. |A^H*A|^2 is a BCCB matrix and admits a fast matrix-vector multiplication
function g = Grad_f(x, S, v, N)
% g = LAH_H( LAH(s, H) - y, H, N );
B = ifft2(S .* fft2( reshape(x, sqrt(N), sqrt(N)) ));
% g = real(B(:)) - v;
g = (B(:) - v);

%% Calculate the gradient of g at the point D^*x_{k-1}
function g = Grad_g(x, D, Dt, Lambda_mu, mu)
y = Dt(x);
g = D(y - Soft(y, Lambda_mu))/mu;

%% Left AH: Implementation of A*Y = H*kron(F, F)*Y efficiently, using FFT operations
function X = LAH(Y, H)
% Determine dimensions
[My, Ny] = size(Y);
[Mh, Nh] = size(H);

X = H * vec(pfft2(reshape(Y, sqrt(My), sqrt(My)), sqrt(Mh), 'fft'));
% X = H * cell2mat( arrayfun(@(ii) vec(pfft2(reshape(full(Y(:, ii)), sqrt(My), sqrt(My)), sqrt(Mh), 'fft')), 1:Ny, 'UniformOutput', false) );

%% Left AH Hermitian: Implementation of A^H*Y = kron(F, F)^H*H^H*Y efficiently, using FFT operations
function X = LAH_H(Y, H, N)
% Determine dimensions
[My, Ny] = size(Y);

Z = H' * Y;

X = cell2mat( arrayfun(@(ii) vec(pfft2_ABZ(reshape(Z(:, ii), sqrt(My), sqrt(My)), sqrt(N), 'fft_h')), 1:Ny, 'UniformOutput', false) );

%% Similar to LAH_H, only the output is a vector
function X = LAH_I(Y, N)
% Determine dimensions
[My, Ny] = size(Y);

X = arrayfun(@(ii) FirstElement(pfft2_ABZ(reshape(full(Y(:, ii)), sqrt(My), sqrt(My)), sqrt(N), 'fft_h')), 1:Ny);

%% Take first element of a matrix
function a = FirstElement(Q)
a = Q(1, 1);

%% Vectorize a matrix
function v = vec( x )
v = x(:);




