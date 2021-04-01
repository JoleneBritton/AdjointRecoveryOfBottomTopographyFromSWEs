%% Accuracy Test 
lr = 0.6; lr_type = 0;
gammaL = 1.1e-6; gammaH = 5e-8; g_hat = 1;
noise_meas = 0; noise_guess = 0; 
k_meas = 3; k_forward = 2; k_adjoint = 1;
final_time = 0.01;
num_of_iters = 1; plot_iters = 100;
accuracy_test = 1;

% Store parameters
hyparams = [lr, lr_type, gammaL, gammaH, final_time, num_of_iters, plot_iters, noise_meas, noise_guess, g_hat];
k = [k_meas, k_forward, k_adjoint];
func_val = [2, 5, 3];
choice = str2double(sprintf('%d',func_val));

% Mesh options
first_iter = 3;
last_iter = 5;
a = 25;
ex_meas_all = (8*a)*2.^[first_iter:last_iter];
ex_forward_all = (a)*2.^[first_iter:last_iter];
ex_adjoint_all = a*2.^[first_iter:last_iter];
nt_meas_all = 100*2.^[first_iter:last_iter];

% Initialize matrices for storing data 
HH = zeros(length(ex_forward_all), k_forward+1, max(ex_forward_all));
QQ = zeros(length(ex_forward_all), k_forward+1, max(ex_forward_all));

for i = 1:length(ex_meas_all)
    % Print info
    disp(['Iteration ', num2str(i), ' of ', num2str(length(ex_forward_all))]);
    
    % Extract mesh data for one run
    ex_forward = ex_forward_all(i);
    ex_adjoint = ex_forward_all(i);
    ex_meas = ex_meas_all(i);
    meshes = [ex_meas, ex_forward, ex_adjoint];
    
    % Num. of time steps
    nt_meas = nt_meas_all(i);
    nt_forward = nt_meas * (ex_forward / ex_meas);
    nt_adjoint = nt_forward * (ex_adjoint / ex_forward);
    nt = [nt_meas, nt_forward, nt_adjoint];
    
    % Load the data class
    data = SWE_Data.data_params(k,nt,func_val,meshes,hyparams,choice,accuracy_test);
    if ~exist(data.front_path, 'dir')
        mkdir(data.front_path)
    end
    
    % Run the forward scheme - if data not yet generated
    data = SWE_Forward(1,data);

    % Load in and store the results
    fileH_T = sprintf('%s/SWE_PredData_%s_H_FinalTime.csv',data.front_path,data.filename_str);
    fileQ_T = sprintf('%s/SWE_PredData_%s_Q_FinalTime.csv',data.front_path,data.filename_str);
    HH(i,:,1:ex_forward) = transpose(csvread(fileH_T));
    QQ(i,:,1:ex_forward) = transpose(csvread(fileQ_T));
end

% Display the errors and order of accuracy
disp('FORWARD DATA')
dx_forward_all = (data.end_pt - data.start_pt)./ ex_forward_all;
SWE_PostProcess.errors(ex_forward_all,HH,QQ,dx_forward_all,final_time,k_forward+1,choice)



