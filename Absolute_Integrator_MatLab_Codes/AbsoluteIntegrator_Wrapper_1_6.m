%                                  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM    
%  Welcome to Absolute Integrator   NMMMMMMMM.                   '=MMMML   
%           Version 1.6.4             MMMMMMMMM                     'MMM   
%           © Lewys Jones              MMMMMMMMM,                    .MM.  
%                                        MMMMMMMMM.                   IM.  
%       Updated 27th May   2014           MMMMMMMMM.                  .M.  
%                                          'MMMMMMMMM.                  '  
%                                           'MMMMMMMMM.                     
%                                             'MMMMMMMMM                   
%             Contact:                          'MMMMMMMM.                  
%   lewys.jones@materials.ox.ac.uk                .MMMMMMMM.               
%                                                   ?MMMMM'                
%              Manual:                             .:MMM'                  
%        www.lewysjones.com                       .MMMM'                   
%                                               .MMMM'                     
%   Requires:                                  .MMMZ'                      
%    - Image Processing Toolbox               MMMM'                      ? 
%    - Curve Fitting Toolbox               .=MMM'                       :M 
%                                         .MMMM'                       MM.
%                                       .MMMM'                        OMMM.
%                                     .cMMM'                      .NMMMMM~ 
%                                   .MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM  
% Workspace setup:                 .=MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM'  
close all, clear all
clc %================= BEGIN User Defined Variables =======================

% Main Setup Options:
loadSimulation              = 0         ; % Load data as if simulation? 1 = yes, 0 = no. If processing data as a simulation, no detector or flux options (next section) apply.

% Detector & Flux Analysis Options:
faradayCurrentImaging       = 1         ; % What was the imaging probe current? (units pA)
faradayCurrentDetector      = 1         ; % What was the detector mapping probe current? (units pA)
weightDetector              = 1         ; % Use 2D detector flux-weighting?  1 = yes, 0 = no. See reference: ADD REFERENCE HERE!!.
    innerAngle              = 45      ; % What is the detector inner angle? (units of milli-radians)
    TDSexponent             = 2.268     ; % If no flux image is available, what exponent should be used? See reference: ADD REFERENCE HERE!!.
analyse_CCD_flux            = 0         ; % Analyse CCD flux image to determine correct flux weighting? 1 = yes, 0 = no.
    convergence_angle       = 21        ; % What is the probe convergence / semi-angle? (units of milli-radians)

% Image Peak-finding Variables:
manual_best_size            = 1         ; % Use a manual overide for feature separation? 0 = automaic, 1 = manual over-ride (specify below).
    best_size               = 19        ; % Size of the manual overide, units of pixels feature separation (must be ODD).
refinePositions             = 1         ; % Refine coordinates to sub-pixel precision? 1 = yes, 0 = no.
removeZeroSpot              = 0         ; % Perform high-pass filter before peak-finding? 0 = no, otherwise units of percentage radial limit to mask FT zero spot.

% Experiment / Sample Specific Variables:
manualPixelWidth            = 0    ; % Zero for automatic (uses diffractogram calibration), otherwise units of nm.
    latticeParameter        = 0.408     ; % Units of nm. Au = 0.408, Pt = 0.392
    calibrationPlanes        = [ 1 1 1 ] ; % Enter the three index spacing notation for the diffractogram calibration e.g. "[ 1 1 1 ]".
subtractBackground          = 0         ; % Subtract an additional background from image, e.g. a carbon support?
    fitLocalBackground      = 1         ; % Fit the background locally? 1 = yes, 0 = use a constant value.
        backfillIterations  = 2000      ; % The number of iterative refinements to use if fitting alocal background.
    add_extra_mask          = 0         ; % Add an extra masking area? 1 = yes, 0 = no.

% Cross-section Integration Options:
upsamplingFactor            = 1         ; % Upsample image to create smoother Voroinoi polgons? 1 = no change, otherwise specify the upsampling factor (>1).
autoRadius                  = 1         ; % Automatically determine the integration cut-off radius? 1 = automatic, zero to over-ride (specify below).
    manualRadius            = 13        ; % Manual integration radius to use (positive real numbers).
    
% Simulation Specific Variables:
simulationTiling            = 1         ; % The number of times the image is tiled in x and y.

% Progress / Display Variables:
show_progress               = 1         ; % Choose whether to display figures during processing.
display_results             = 1         ; % Choose whether to display figures after processing.

%======================== BETA / DEBUG FEATURES ===========================
use_SpectrumIntegrator      = 0         ; % NEW feature - testing only.
analyseAzimuthalSensitivity = 0         ;
calculateEllipticity        = 0         ; % WARNING! May take > 1 hour for large images.
offsetDetector              = [ 0 0 ]   ; % [Y X] Percetage of inner angle to offset the detector.
addRandomOffset             = 0         ; % Adds a random zero-mean offset to the peak positions with this standard deviation.
addNoiseTargetSNR           = 0         ; % Zero for noise-free, otherwise positive.
%====================== END User Defined Variables ========================
tic
version_string = ('AbsoluteIntegrator v1.6.4') ;
run AbsoluteIntegrator_Core_1_6_4
toc