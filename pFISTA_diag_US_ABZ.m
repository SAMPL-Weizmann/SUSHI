function [ x ] = pFISTA_diag_US_ABZ( y, H, Params )
% PFISTA_DIAG_US - Performs the Fast Proximal Gradient method on the problem
% 
%              min_{x} \lambda|| x ||_1 + 0.5|| y - Ax ||_2^2
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
% [ x ] = pFISTA_diag_US_2( y, H, Params )
%
% Inputs:
% -------
% y      - Input M^2 X 1 observations matrix (Generaly complex)
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
% Written by Oren Solomon, Technion I.I.T. Ver 1. 27-09-2016
% Revised by Avinoam Bar-Zion, Caltech Ver 2. 02-10-2017

global VERBOSE

%% Flags
Calc_S_flag = 0; 

%% Initializations
% -----------------------------------------------------------------------------------------------------------------------
% Paramters
beta       = Params.Beta;
L          = Params.L0;              % Lipschitz constant
lambda_bar = Params.LambdaBar;
lambda     = Params.Lambda;          % l1 regularization parameters
N          = Params.N;               % Length of x

% Init
t          = 1;
t_prev     = 1;

% Memory allocation
x          = zeros(N, 1);
x_prev     = x;

% Vectorize
y = y(:);

if VERBOSE; disp('pFISTA_diag: Calculating the Lipschitz constant...'); end;

% Calculate S 
if VERBOSE; fprintf('S_calc: ');t2 = tic; end;
if Calc_S_flag
    S = S_calc(H, N); % ???
else
    load('Smat_64_to_512.mat');
end
if VERBOSE; toc(t2); end;

if VERBOSE; fprintf('v_calc: ');t3 = tic; end;
v = v_calc(y, H, N); % ???
if VERBOSE; toc(t3); end;

if isempty(L) % 2*Lambda_max(A^H*A) according to FISTA paper
    L = 2*real(max(max(S))); % ??? is it the same as in the fista paper?
end

if VERBOSE; disp('done.'); end;

%% Iterations
% -----------------------------------------------------------------------------------------------------------------------
if VERBOSE; disp('pFISTA_diag: Running iterations...'); end;
for kk = 1:Params.IterMax
    Titer = tic;
    if VERBOSE; fprintf(['Iteration #' num2str(kk) ': ']); end;
    
    % Z update
    z = x + ((t_prev - 1)/t)*(x - x_prev);  
    
    % Gradient step: G = Z - (1/L)*( A^H*(A*x - y) ) % ??? I am not sure about
    % the (1/L) should read again FISTA
    g = z - (2/L)*Grad(z, S, v, N);

    % Soft thresholding
    x_prev = x;
    x = sign(g).*max(abs(g) - lambda/L, 0);
    
    % Projection onto the non-negative orthant, only for a non-negative constraint
    if Params.NonNegOrth == 1 % should work after fixing the bug
%         x = real(x);
%         x(x < 0) = 0;
%         x = (x);
        x(abs(x) < 0) = 0;
    end
    
    
    % Parameter updates for next iteration
    t_prev = t;
    t = 0.5*(1 + sqrt(4*t^2 + 1));
    
    lambda = max(beta*lambda, lambda_bar); % From Yonina's matrix FISTA paper. do we need this ???
    
    if VERBOSE; toc(Titer); end;
end
if VERBOSE; disp('Done pFISTA_diag.'); end;

function v = v_calc(y, H, N) % ??? what is the difference between hermit and transpose
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
function g = Grad(x, S, v, N)
% g = LAH_H( LAH(x, H) - v, H, N );
B = ifft2(S .* fft2( reshape(x, sqrt(N), sqrt(N)) ));

g = (B(:) - v);

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
