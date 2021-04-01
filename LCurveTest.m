%% L-curve Test- using previously generated data
lr = 0.6; lr_type = 0;
gammaL = 1e-6; gammaH = 5e-8; 
noise_meas = 0.1; noise_guess = 0.25;
ex_meas = 400; ex_forward = 50; ex_adjoint = 25;
k_meas = 3; k_forward = 2; k_adjoint = 1;
nt_meas = 5000;
nt_forward = nt_meas * (ex_forward / ex_meas);
nt_adjoint = nt_forward * (ex_adjoint / ex_forward);
final_time = 0.05;
num_of_iters = 1000; plot_iters = 100;
accuracy_test = 0;

% Store parameters
meshes = [ex_meas, ex_forward, ex_adjoint];
k = [k_meas, k_forward, k_adjoint];
nt = [nt_meas, nt_forward, nt_adjoint];
func_val = [2, 5, 3];
g_hat_all = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000];

% Initialize vectors for storing results
err = zeros(length(g_hat_all),1);
reg = zeros(length(g_hat_all),1);

i = 1;
for g_hat = g_hat_all
    disp(g_hat);
    hyparams = [lr, lr_type, gammaL, gammaH, final_time, num_of_iters, plot_iters, noise_meas, noise_guess, g_hat];
    
    choice = str2double(sprintf('%d',func_val));
    data = SWE_Data.data_params(k,nt,func_val,meshes,hyparams,choice,accuracy_test);
    filename_str = data.filename_str;
    
    % Read in and store data
    file_errors = sprintf('%s/SWE_IterationErrors_%s.csv',data.front_path,data.filename_str);
    file_besterror = sprintf('%s/SWE_BestError_%s.csv',data.front_path,data.filename_str);
    errors = csvread(file_errors);
    best_err = csvread(file_besterror);
    err(i) = errors(best_err);
    
    file_bestreg = sprintf('%s/SWE_BestReg_%s.csv',data.front_path,data.filename_str);
    reg(i) = csvread(file_bestreg);
        
    i = i+1;
end

% Use Latex formatting
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultTextInterpreter', 'latex');

% Plotting Line Styles
lstyles = ["--", ":", "-.", "-"];

% Font & Marker Sizes
title_fs = 85;
subtitle_fs = 50;
xlabel_fs = 75;
legend_fs = 75;
sublegend_fs = 35;
ticks_fs = 50;
marker_s1 = 800;
marker_s2 = 250;

% Plot
figure;
loglog(err, reg, '-o', 'LineWidth', 3, 'MarkerSize', 20)
xlabel('Residual Error', 'FontSize', xlabel_fs);
ylabel('Magnitude of Regularizer', 'FontSize', xlabel_fs);
ax = gca;
ax.XAxis.FontSize = ticks_fs;
ax.YAxis.FontSize = ticks_fs; 
set(gcf,'units','normalized','outerposition',[0 0 0.75 1])
    
