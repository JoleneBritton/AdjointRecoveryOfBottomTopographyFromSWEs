function SWE_Main(k,nt,func_val,meshes,hyparams,accuracy_test)

    %% ---------------------------------------------------------------------------------------
    %     Set Up
    % ----------------------------------------------------------------------------------------

    close all;
    clc;
    
    % Start timer
    tic;
    
    % Get info
    choice = str2double(sprintf('%d',func_val));
    
    % Initialize and fill data class
    data = SWE_Data.data_params(k,nt,func_val,meshes,hyparams,choice,accuracy_test);
    if ~exist(data.front_path, 'dir')
        mkdir(data.front_path)
    end
    
    % Display info
    display(choice);
    display(hyparams);
    display(meshes);
    display(data.filename_str);
    display(data.p_true);
    display(data.p_guess);
        
    % Generate measured data
    data.meas_data = 1;
    data = SWE_Meas(data);

    % Initialize iteration counter
    m_iter = 1;
    data.best_err1 = 1000000000;
    data.meas_data = 0;

    % Iterate
    while m_iter <= data.num_of_iters
        % Forward problem
        data = SWE_Forward(m_iter,data);

        % Initial choice for Z is initial guess for P
        if m_iter == 1
            Z = data.p_guess;
        end
        
        % Inverse problem to update bottom function
        [data, Z, m_iter] = SWE_Inverse.inverse(Z,m_iter,data);
    end
            
    % Save data
    file_errors = sprintf('%s/SWE_IterationErrors_%s.csv',data.front_path,data.filename_str);
    file_besterror = sprintf('%s/SWE_BestError_%s.csv',data.front_path,data.filename_str);
    file_bestreg = sprintf('%s/SWE_BestReg_%s.csv',data.front_path,data.filename_str);
    file_errorsreg = sprintf('%s/SWE_IterationErrorsReg_%s.csv',data.front_path,data.filename_str);
    file_bestp = sprintf('%s/SWE_BestP_%s.csv',data.front_path,data.filename_str);
    file_bestH = sprintf('%s/SWE_BestH_%s.csv',data.front_path,data.filename_str);
    file_bestQ = sprintf('%s/SWE_BestQ_%s.csv',data.front_path,data.filename_str);
    file_bestB = sprintf('%s/SWE_BestB_%s.csv',data.front_path,data.filename_str);
    
    dlmwrite(file_errors, data.iter_err1,'delimiter',',','precision',20)
    dlmwrite(file_besterror, data.best_iter1,'delimiter',',','precision',20)
    dlmwrite(file_bestreg, data.reg,'delimiter',',','precision',20)
    dlmwrite(file_errorsreg, data.iter_err_reg,'delimiter',',','precision',20)
    dlmwrite(file_bestp, data.best_p1,'delimiter',',','precision',20)
    dlmwrite(file_bestH, data.H_all_pred_best1,'delimiter',',','precision',20)
    dlmwrite(file_bestQ, data.Q_all_pred_best1,'delimiter',',','precision',20)
    dlmwrite(file_bestB, data.B_all_pred_best1,'delimiter',',','precision',20)
    
    % Print basic stats
    disp(sprintf('Best Iteration: %d', data.best_iter1))
    disp(sprintf('Best Error:     %d', data.best_err1))


end
    
    
        