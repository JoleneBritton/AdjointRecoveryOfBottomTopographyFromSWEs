% n = case number you want to run

% run_type = 0; % run main
% run_type = 1; % run results
% run_type = 2; % run both

function [] = SWE_Run(n, run_type)
n=floor(n);
display(n);

%% ---------------------------------------------------------------------------------------
%     Global Parameters
% ----------------------------------------------------------------------------------------

% Turn off accuracy test_settings
accuracy_test = 0;

% Learning rate and type
lr = 0.6; lr_type = 0;


%% ---------------------------------------------------------------------------------------
%     Test Dependent Parameters
% ----------------------------------------------------------------------------------------

switch n
    case 1
        %% CASE 271/4.4b/Fig 16: ptrue = 2 bumps and 1 well, pguess = constant 1
        % Regularization parameters
        gammaL = 1e-5; gammaH = 1e-6; g_hat = 1;
        
        % Noise sizes
        noise_meas = 0.1; noise_guess = 0.25; 
        
        % Mesh sizes - spatial
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25; 
        
        % Temporal parameters
        nt_meas = 10000; final_time = 0.2;
        
        % Piecewise polynomial degrees
        k_meas = 3; k_forward = 2; k_adjoint = 1;

        % Iteration parameters
        num_of_iters = 1000; plot_iters = 100;

        % Determines ICs, p_true, p_guess
        func_val = [2, 7, 1];

    case 2
        %% CASE 211/4.4a/Fig 14: ptrue = one centered bump, pguess = constant 1
        gammaL = 1e-5; gammaH = 1e-6; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25; 
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 10000; final_time = 0.2;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 1, 1];

    case 3
        %% CASE 273/4.1f/Fig 6: ptrue = 2 uneven bumps (flipped), pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-9; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25; 
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 7, 3];

    case 4
        %% CASE 263/4.1e/Fig 5: ptrue = 2 uneven bumps (flipped), pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 6, 3];

    case 5
        %% CASE 258/4.2b/Fig 9: ptrue = 2 uneven bumps (flipped), pguess = 1 bump function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.2; %0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25; 
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 10;
        func_val = [2, 5, 8];

    case 6
        %% CASE 257/4.2c/ Fig 10: ptrue = 2 uneven bumps, pguess = well
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25; 
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 7];

    case 7
        %% CASE 253/4.3c/Fig 13: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 0; gammaH = 0; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000;final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 8
        %% CASE 253/4.3b/Fig 12: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 0; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 5000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 9
        %% CASE 253/4.3a/Fig 11: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 0; gammaH = 5e-8; g_hat = 1; 
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 10
        %% CASE 253/L-curve1: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1e-5;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 11
        %% CASE 253/L-curve2: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1e-4;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 12
        %% CASE 253/L-curve3: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1e-3;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 13
        %% CASE 253/L-curve4: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1e-2;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 14
        %% CASE 253/L-curve5: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1e-1;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 15
        %% CASE 253/L-curve6: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 10;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 16
        %% CASE 253/L-curve7: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 100;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 17
        %% CASE 253/L-curve8: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1000;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 18
        %% CASE 253/L-curve9: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 10000;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 19
        %% CASE 253/L-curve10: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 100000;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 20
        %% CASE 253/Mesh1: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 25; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 2;
        nt_meas = 5000;  final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 21
        %% CASE 253/Mesh2: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 50;
        k_meas = 3; k_forward = 2; k_adjoint = 2;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 22
        %% CASE 253/4.1d(b)/Fig 4(b): ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 23
        %% CASE 253/4.1d(a)/Fig 4(a): ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 1e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 3];

    case 24
        %% CASE 251/4.2a/Fig 8: ptrue = 2 uneven bumps (flipped), pguess = constant function
        gammaL = 1e-6; gammaH = 1e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 5, 1];

    case 25
        %% CASE 244/4.1c/Fig 3: ptrue = 2 even bumps, pguess = 1 tall centered bump
        gammaL = 1e-6; gammaH = 1e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 4, 4];

    case 26
        %% CASE 235/4.1b/Fig 2: ptrue = right bump, pguess = left bump
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 3, 5];

    case 27
        %% CASE 226/4.1a/Fig 1: ptrue = left bump, pguess = right bump
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25; 
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25; 
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1000; plot_iters = 100;
        func_val = [2, 2, 6];
        
    case 28
        %% CASE 253/ConvergencePlotTest: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 1e-6; gammaH = 5e-8; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1500; plot_iters = 300;
        func_val = [2, 5, 3];
        
    case 29
        %% CASE 253/ConvergencePlotTest: ptrue = 2 uneven bumps, pguess = oscillatory function
        gammaL = 5e-6; gammaH = 1e-7; g_hat = 1;
        noise_meas = 0.1; noise_guess = 0.25;
        ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
        k_meas = 3; k_forward = 2; k_adjoint = 1;
        nt_meas = 5000; final_time = 0.05;
        num_of_iters = 1500; plot_iters = 300;
        func_val = [2, 5, 3];
end


%% ---------------------------------------------------------------------------------------
%     Run &/or Plot
% ----------------------------------------------------------------------------------------

% Num. of timesteps 
nt_forward = nt_meas * (ex_forward / ex_meas);
nt_adjoint = nt_forward * (ex_adjoint / ex_forward);
        
% Store parameters
hyparams = [lr, lr_type, gammaL, gammaH, final_time, num_of_iters, plot_iters, noise_meas, noise_guess, g_hat];
meshes = [ex_meas, ex_forward, ex_adjoint];
k = [k_meas, k_forward, k_adjoint];
nt = [nt_meas, nt_forward, nt_adjoint];

% Run main code
if (run_type == 0) || (run_type == 2)
    SWE_Main(k,nt,func_val,meshes,hyparams,accuracy_test)
end
close all

% Generate plots
if (run_type == 1) || (run_type == 2)
    SWE_Results(k,nt,func_val,meshes,hyparams,accuracy_test)
end
end
