% ----------------------------------------------------------------- %
% This script demonstrates sparse deconvolution on SOFI processed
% images using the FISTA algorithm. The FISTA codes are taken from
% Amir beck's freely available codes.
% ----------------------------------------------------------------- %

clc;
clear;
% close all;

global VERBOSE
VERBOSE = 1;

% Display images to screen
DisplayFlag = 1;

if VERBOSE; disp('------------------------------------------------------'); end;
if VERBOSE; disp('~                    SUSHI Ver. 1                    ~'); end;
if VERBOSE; disp('------------------------------------------------------'); end;
if VERBOSE; disp(' '); end;

% -------------------------------------------------------------------------
%% Add relevant folders and subfolders
% -------------------------------------------------------------------------
% addpath(genpath('D:\SUSHI\150Frames_dual'));
addpath(genpath('E:\OREN_DATA\CS-DCEUS\Datasets\Raw_data'));
addpath(genpath('D:\MATLAB\Papers\Paper_13_Nat_SUSHI\rwt-master'));
addpath(genpath('D:\SUSHI\SUSHI_Code\pFISTA_diag'));
addpath(genpath('..\Patch data\Data'));               % Raw data
addpath(genpath('..\Patch data\HNO'));                % Utility functions
addpath(genpath('..\Patch data\FISTA'));              % FISTA
addpath(genpath('..\Patch data\FISTA_TV'));           % TV FISTA
addpath(genpath('..\Patch data\rwt-master'));         % Rice wavelet toolbox V.3
addpath(genpath('..\Patch data\iLET'));               % iLET 
addpath(genpath('..\Patch data\pFISTA_diag'));        % CS FISTA 

% -------------------------------------------------------------------------
%% load raw data for processing
% --------------------------------------------------------------------------
DataTypeFlag = 'mat';
PatchNum     = 7;
switch PatchNum
    case 1
        DataTypeFolder = 'Raw_data\Patch1';
        DataTypName    = 'Patch1_128X128_LowDens';
    case 2
        DataTypeFolder = 'Raw_data\Patch2';
        DataTypName    = 'Patch2_128X128';
    case 3
        DataTypeFolder = 'Raw_data\Patch3';
        DataTypName    = 'Patch3_128X128_HighDens';
    case 4
        DataTypeFolder = 'Raw_data\Patch4';
        DataTypName    = 'Patch4_128X128_HighDens';
    case 5
        DataTypeFolder = 'Raw_data\Patch5';
        DataTypName    = 'Patch5_128X128_HighDens';
    case 6
        DataTypeFolder = 'Raw_data\Patch6';
        DataTypName    = 'Patch6_calibration';
    case 7 % best patch
        DataTypeFolder = 'Raw_data\Patch7';
        DataTypName    = 'Patch7_calibration';
    case 8
        DataTypeFolder = 'Raw_data\Patch8';
        DataTypName    = 'Patch8_128X128';
    case 9
        DataTypeFolder = 'Raw_data\Patch9';
        DataTypName    = 'Patch9_64X64';
    case 10
        DataTypeFolder = 'Raw_data\Patch10';
        DataTypName    = 'Patch10_64X64';
    case 11
        DataTypeFolder = 'Raw_data\Patch11';
        DataTypName    = 'Patch11_64X64';
    case 12
        DataTypeFolder = 'Raw_data\Patch12';
        DataTypName    = 'Patch12_64X64';
    case 13
        DataTypeFolder = 'Raw_data\Patch13';
        DataTypName    = 'Patch13_64X64';
    case 14
        DataTypeFolder = 'Raw_data\Patch14';
        DataTypName    = 'Patch14_64X64';
    case 15
        DataTypeFolder = 'Raw_data\Patch15';
        DataTypName    = 'Patch15_64X64';
    case 16
        DataTypeFolder = 'Raw_data\Patch16';
        DataTypName = 'Patch16_64X64';
    case 17
        DataTypeFolder = 'Raw_data\Patch17';
        DataTypName = 'Patch17_64X64';
    case 18
        DataTypeFolder = 'Raw_data\Patch17';
        DataTypName = 'BrainPatchMid'; 
    case 19
        DataTypeFolder = 'Raw_data\Patch17';
        DataTypName = 'BrainPatchCortex';
    case 20
        DataTypeFolder = 'Raw_data\Patch17';
        DataTypName = 'BrainPatchFunction';  
    case 21
        DataTypeFolder = 'Raw_data\Patch17';
        DataTypName = 'BrainPatchCortex2'; 
    case 22
        DataTypeFolder = 'Raw_data\Patch17';
        DataTypName = 'LongKidney';
    case 23
        DataTypeFolder = 'Raw_data\Patch17';
        DataTypName = 'fUS';
end


% Load Movie
if VERBOSE; fprintf('Loading data...'); end;
switch lower(DataTypeFlag)
    case 'mat'
        % Assume that the new data structure is called 'Patch' and it is a
        % struct with several fields
%         load(fullfile(DataTypeFolder, DataTypName));
        load(fullfile(DataTypName));
        frameCnt = 5;
        Xs    = 1;
        Ys    = 1;
        Dsize = 64;
%         Patch.AC2_pos      = Patch.AC2_pos(:,:,frameCnt);
%         Patch.AC2_neg      = Patch.AC2_neg(:,:,frameCnt);
%         Patch.AC2_combined = Patch.AC2_combined(:,:,frameCnt);
        
        % Determine dimensions
        [M, N, K] = size(Patch.AC2_neg);
    otherwise
        error('Unknown data type.');
end
if VERBOSE; disp('Done.'); end;

% -------------------------------------------------------------------------
%% Perform sparse deconvolution on data
% -------------------------------------------------------------------------
% Obsolete
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Deconvolution type
% FISTAType = 'l1';
% FISTAType = 'wavelet';
% FISTAType = 'dct';
% FISTAType = 'tv';
% FISTAType = 'ilet';
% FISTAType = 'lucy';
%
% % Parameters
% lambda         = 0.0001;           % Regularization: try ~300 for l1 ~1 for TV and ~1 for wavelet
% params.MAXITER = 100;              % Maximum number of iterations
% params.fig     = 0;                % 0 - do not display images from within the FISTA functions
% params.BC      = 'periodic';       % Boundary conditions - 'periodic' or 'reflexive'
% 
% % For TV deblurring
% params.denoiseiter = 50;           % Number of iterations per denoising loop (prox estimation)
% params.TV          = 'iso';        % 'iso' for isotropic or 'l1' for anisotropic TV minimization
% LowBound           = 0;            % Lower bound constrain for X_out
% UpBound            = Inf;          % Upper bound constrain for X_out
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Choose solver
% -------------------------------------------------------------------------
FISTAType = 'pFISTA_diag';
% FISTAType = 'pFISTA_analysis';
% FISTAType = 'pFISTA_diag_TV';


% Super-resolution factor
% -------------------------------------------------------------------------
SRF                  = 8;          

% Solver parameters
% -------------------------------------------------------------------------
% Parameters to choose (SRF = 8):
% l1     : lambda = 0.5, IterMax = 500; 0.12/500 for Patch 15
% wave2d : lambda = 0.2, IterMax = 100
% TV     : Does not seem to perform well

% For SR using pFISTA_diag_US
cs_params.Beta       = 1;
cs_params.L0         = [];         % Lipschitz constant
cs_params.LambdaBar  = 1e-7;
cs_params.Lambda     = 1; 0.5; 0.5; 1;        % 1 - l1 regularization parameters. 0.0001 for TV, 0.009 for analysis
cs_params.N          = (N*SRF)^2;  % Length of x
cs_params.IterMax    = 500;
cs_params.NonNegOrth = 1;
cs_params.MaxTimeLag = 1;

% Additional for pFISTA_TV_US
cs_params.DenoiseIter = 50;
cs_params.TV          = 'iso';
cs_params.mon         = 0;

% Additional for pFISTA_analysis_US
cs_params.AnalysisType = 'wave2d';
cs_params.mu           = 1e-4;

% Run Solver
% -------------------------------------------------------------------------
if VERBOSE; disp(['Runnig solver: ' FISTAType]); end;
if VERBOSE; disp(['-------------------------------']); end;
CurrentMovie = {'neg', 'pos', 'combined'};
for kk = 1:length(CurrentMovie) 
    if VERBOSE; disp(' '); end;
    if VERBOSE; disp(['Filename: ' DataTypName ', Direction: ' CurrentMovie{kk} ' flow.']); end;
    if VERBOSE; disp(['~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~']); disp(' '); end;
    
    % Choose datafile to work on
    ChoiceFlag = CurrentMovie{kk};
    switch lower(ChoiceFlag)
        case 'pos'      % Positive flow
            CurrPatch        = Patch.pos;
            AC2              = Patch.AC2_pos;
%             cs_params.Lambda = 0.5;
        case 'neg'      % negative flow
            CurrPatch        = Patch.neg;
            AC2              = Patch.AC2_neg;
%             cs_params.Lambda = 0.5;
        case 'combined' % Unfiltered
            CurrPatch        = Patch.dl;
            AC2              = Patch.AC2_combined;
%             cs_params.Lambda = 0.5;
    end
    
    % PSF^2
    psfSquared = Patch.PSF2;
    
    % Choose solver
    switch lower(FISTAType)
        % "Regular" l1 sparsity based deconvolution
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case 'l1'
            D = @(x) x;
            [X_out, fun_all] = deblur_wavelet_FISTA_trans(AC2, psfSquared, center, D, D, lambda, params);
            
            % Wavelet sparsity based deconvolution
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case 'wavelet'
            WaveD     = 2;                                        % Wavelet depth hierarchy - Daubechies wavelets
            WaveScale = daubcqf(32, 'mid');                       % Type of wavelet filter
            D         = @(x) midwt(x, WaveScale, WaveD);          % Synthesis operator
            Dt        = @(x) mdwt(x, WaveScale, WaveD);           % Analysis operator
            [X_out, fun_all] = deblur_wavelet_FISTA_trans(AC2, psfSquared, center, D, Dt, lambda, params);
            
            % DCT sparsity based deconvolution
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case 'dct'
            D  = @(x) dct2(x);
            Dt = @(x) idct2(x);
            [X_out, fun_all] = deblur_wavelet_FISTA_trans(AC2, psfSquared, center, D, Dt, lambda, params);
            
            % Total variation sparsity based deconvolution
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case 'tv'
            [X_out, fun_all] = deblur_tv_fista(AC2, psfSquared, center, lambda, LowBound, UpBound, params);
            
            %iLET sparsity based deconvolution
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case 'ilet'
            [X_out, ~, fun_all, ~, ~, ~, ~] = iLET_deconv(AC2, psfSquared, 'showvalue', 1, 'LAMBDA', lambda);  % Doesn't seem to work
            
            % Richardson-Lucy deconvolution
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case 'lucy'
            X_out = deconvlucy(AC2, psfSquared);
            fun_all = [];
            
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Super-resolution via compressed sensing
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case 'pfista_diag'
            % Convert data and PSF to the frequency domain
            w = tukeywin(64,0.5);
            for ff = 1:K
                
                AC_F(:, :, ff) = ((1/N^2)*fft2((AC2(:, :, ff).*(w*w'))));
            end
            psfSquared_F = ((1/N^2)*fft2((psfSquared)));
            
            % Perform column normalization on A
            Nfactor = norm(psfSquared_F(:), 2); % ???
            
            % Run algorithm
%             X_out = pFISTA_diag_US_3( reshape(AC_F, M*N, K), diag(psfSquared_F(:)/Nfactor), cs_params );
            AC_F = AC_F(:, :, 1); X_out = pFISTA_diag_US_ABZ( AC_F(:), diag(psfSquared_F(:)/Nfactor), cs_params );
%                     X_out = pFISTA_diag_US( AC_F(:), diag(psfSquared_F(:)/Nfactor), cs_params );
            fun_all = [];
            
            % Output as an image
            X_out = (real(reshape(X_out, sqrt(cs_params.N), sqrt(cs_params.N), cs_params.MaxTimeLag)));
        case 'pfista_diag_tv'
            % Convert data and PSF to the frequency domain
            for ff = 1:K
                AC_F(:, :, ff) = fftshift((1/N^2)*fft2(AC2(:, :, ff)));
            end
            psfSquared_F = fftshift((1/N^2)*fft2(psfSquared));
            
            % Perform column normalization on A
            Nfactor = norm(psfSquared_F(:), 2);
            
            % Run algorithm
            AC_F = AC_F(:, :, 1); X_out = pFISTA_diag_TV_US( AC_F(:), diag(psfSquared_F(:)/Nfactor), cs_params );
            fun_all = [];
            
            % Output as an image
            X_out = fftshift(real(reshape(X_out, sqrt(cs_params.N), sqrt(cs_params.N))));
        case 'pfista_analysis'
            % Convert data and PSF to the frequency domain
            w = tukeywin(64,0.5);
            for ff = 1:K
                
                AC_F(:, :, ff) = ((1/N^2)*fft2((AC2(:, :, ff).*(w*w'))));
            end
            psfSquared_F = ((1/N^2)*fft2(psfSquared));
            
            % Perform column normalization on A
            Nfactor = norm(psfSquared_F(:), 2);
            
            % Run algorithm
            AC_F = AC_F(:, :, 1); X_out = pFISTA_diag_analysis_US_ABZ( AC_F(:), diag(psfSquared_F(:)/Nfactor), cs_params );
            fun_all = [];
            
            % Output as an image
            X_out = (real(reshape(X_out, sqrt(cs_params.N), sqrt(cs_params.N))));
        otherwise
            error('FISTA Type not supported.');
    end
    
    % Rearrange for display
    % ---------------------------------------------------------------------
    switch lower(ChoiceFlag)
        case 'pos'      % Positive flow
            Patch.SUSHI_pos      = abs(X_out);
        case 'neg'      % negative flow
            Patch.SUSHI_neg      = abs(X_out);
        case 'combined' % Unfiltered
            Patch.SUSHI_combined = abs(X_out);
    end
end
% Combine flows
P = 0.7;
Patch.SUSHI_Overlay = imfuse(Patch.SUSHI_pos.^P, Patch.SUSHI_neg.^P, 'falsecolor');
Patch.SUSHI_overlay = Patch.SUSHI_pos - Patch.SUSHI_neg;
if VERBOSE; disp(['----------------------------------']); end;
if VERBOSE; disp(['Finished runnig solver: ' FISTAType]); end;

% -------------------------------------------------------------------------
%% Display images
% -------------------------------------------------------------------------
if DisplayFlag
    P = 0.4;
    figure; colormap(gray.^P);
    % First row
%     h(1) = subplot(3, 2, 1); 
%     imagesc(imresize(Patch.pos_dl, size(Patch.SUSHI_pos))); axis square; title('Diffraction limited - positive flow');
    h(2) = subplot(3, 2, 1);
    AC2_pos_resize = imresize(Patch.AC2_pos, size(Patch.SUSHI_pos));
    imagesc(AC2_pos_resize(10:end-10, 10:end-10, 1)); axis square; title('AC2');
    h(3) = subplot(3, 2, 2);
    SUSHI_pos_filt = medfilt2(Patch.SUSHI_pos);%imfilter(Patch.SUSHI_pos, psfSquared.^P);
    imagesc(SUSHI_pos_filt(10:end-10, 10:end-10, 1)); axis square; title(['SUSHI, \lambda = ' num2str(cs_params.Lambda)]);
    
    % Second row
%     h(4) = subplot(3, 2, 4);
%     imagesc(imresize(Patch.neg_dl, size(Patch.SUSHI_neg))); axis square; title('Diffraction limited - negative flow');
    h(5) = subplot(3, 2, 3);
    AC2_neg_resize = imresize(Patch.AC2_neg, size(Patch.SUSHI_neg));
    imagesc(AC2_neg_resize(10:end-10, 10:end-10, 1)); axis square; title('AC2');
    h(6) = subplot(3, 2, 4);
    SUSHI_neg_filt = medfilt2(Patch.SUSHI_neg);%imfilter(Patch.SUSHI_neg, psfSquared.^P);
    imagesc(SUSHI_neg_filt(10:end-10, 10:end-10, 1)); axis square; title(['SUSHI, \lambda = ' num2str(cs_params.Lambda)]);
    
    % Third row
%     h(7) = subplot(3, 3, 7);
%     imagesc(imresize(Patch.tot_dl, size(Patch.SUSHI_combined))); axis square; title('Diffraction limited - unfiltered');
    h(8) = subplot(3, 2, 5);
    AC2_combined_resize = imresize(Patch.AC2_combined, size(Patch.SUSHI_combined));
    imagesc(fliplr(AC2_combined_resize(10:end-10, 10:end-10, 1))); axis square; title('AC2');
    h(9) = subplot(3, 2, 6);
    SUSHI_combined_filt = imfilter(Patch.SUSHI_combined, psfSquared.^P);
    imagesc(SUSHI_combined_filt(10:end-10, 10:end-10, 1)); axis square; title(['SUSHI, \lambda = ' num2str(cs_params.Lambda)]);
%     linkaxes(h, 'xy');

    figure    
    h(5) = subplot(2, 1, 1);
    imagesc(AC2_neg_resize(270:456, 235:421, 1)); axis square; title('AC2');
    h(6) = subplot(2, 1, 2);
    imagesc(SUSHI_neg_filt(270:456, 235:421, 1)); axis square; title(['SUSHI, \lambda = ' num2str(cs_params.Lambda)]);
    colormap(gray.^0.6);  
    
    figure
    plot(SUSHI_neg_filt(352, 235:421, 1))
    
    figure
    subplot(2,1,2)
    imshow(cat(3,mat2gray(SUSHI_neg_filt(24:end-24, 24:end-24, 1).^0.6,[0 0.15]), mat2gray(SUSHI_pos_filt(24:end-24, 24:end-24, 1).^0.6,[0 0.15]), mat2gray(SUSHI_pos_filt(24:end-24, 24:end-24, 1).^0.6,[0 0.15])),[])
    axis square; title('Nat SUSHI');
    subplot(2,1,1)
        imagesc(fliplr(AC2_combined_resize(24:end-24, 24:end-24, 1))); axis square; title('Power Doppler');
    colormap(gray.^0.6); 
    axis off
    % Display overlay
    % ---------------------------------------------------------------------
%     figure; colormap(gray.^0.6);
% %     q(1) = subplot(1, 3, 1);
% %     imagesc(imresize(Patch.Overlay, size(Patch.SUSHI_combined))); axis square; title('diffraction limited, combined flows');
%     q(2) = subplot(1, 2, 1);
%     AC2_combined_resize = imresize(Patch.AC2_combined, size(Patch.SUSHI_combined));
%     imagesc(AC2_combined_resize(10:end-10, 10:end-10, 1)); axis square; title('AC2, combined flows');
%     q(3) = subplot(1, 2, 2);
%     SUSHI_Overlay_filt = imfilter(Patch.SUSHI_Overlay, psfSquared.^P);
%     imagesc(SUSHI_Overlay_filt(10:end-10, 10:end-10, :)); axis square; title('SUSHI, combined flows');
%     linkaxes(q, 'xy');
end

%% Save output
% -------------------------------------------------------------------------
% save(fullfile(DataTypeFolder, DataTypName), 'Patch');

% -------------------------------------------------------------------------
%% Save reconstructions
% -------------------------------------------------------------------------
% % Save figures
% I{1} = MeanImage;
% I{2} = abs(dataResizePatch(:, :, 43));
% 
% % STORM
% TS = double(sum(imread('ThunderSTORM\Patch1\hist.tif', 'tif'), 3));
% I{3} = imfilter(imresize(TS, [512 512], 'method', 'bilinear'), psfSquared.^1);
% 
% % SOFI
% I{4} = AC2;
% I{5} = AC4;
% 
% % SUSHI
% I{6} = imfilter(X_out, psfSquared.^1);

% save Images I;
% -------------------------------------------------------------------------

%% filtering signal

% h = fwind2(Hd, win)