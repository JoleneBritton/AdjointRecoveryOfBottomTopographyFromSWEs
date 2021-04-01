function data = SWE_Forward(m_iter, data)

    %% -----------------------------------------------------------------------------------
    %     Parameters & Initial Conditions
    % ------------------------------------------------------------------------------------
    
    % Unpack variables
    ex = data.forward_ex;
    dx = data.forward_dx;
    x_center = data.forward_mesh_center;
    pd = data.forward_k+1;

    data.pd = pd;
    data.ex = ex;
    data.dx = dx;
    data.x_center = x_center;
    bc_type = data.bc_type;
    vander_pd = data.vander_pd;
    
    % Information for P
    p_guess = data.p_guess;

    % Initial conditions for H,Q & Determine B
    [H,Q,B0,B1] = SWE_SetUp.initial_conditions(data,pd,ex,dx,x_center);
    data.B1 = B1;

    et = 500;

    % Initialize matrices for holding all data for space and time
    H_all = zeros(ex,data.pd,et);
    [Q_all, B_all] = deal(H_all);
    
    H_all(:,:,1) = H;
    Q_all(:,:,1) = Q;
 
    
    %% -----------------------------------------------------------------------------------
    %     Begin DG Method
    % ------------------------------------------------------------------------------------

    % Current time
    time = 0;

    % Run DG scheme when final time is greater than 0
    if data.T ~= 0

        % Time step index
        data.fr = 1;  
            
        % Step through time steps
        for timestep = 1:data.forward_nt

            % B.C. padding for initial functions
            Hn = SWE_SetUp.bc_padding(H,bc_type);
            Qn = SWE_SetUp.bc_padding(Q,bc_type);
            B0 = SWE_SetUp.bc_padding(B0,bc_type);
            B1 = SWE_SetUp.bc_padding(B1,bc_type);
            
            % Rename variables
            H0 = Hn;
            Q0 = Qn;

            % 3-Step Runge-Kutta
            for rk_step = 1:3

                %% --------------------------------------------------------
                %     Find dt & Runge-Kutta Steps
                % ---------------------------------------------------------
                    
                if rk_step == 1
                    % Time step size & alpha calc
                    if  m_iter == 1
%                         if data.accuracy_test == 1
                            data.alpha = SWE_SetUp.swe_alpha(H0,Q0,data);
%                         end
                        data.dt = data.forward_dt;
                    else
                        data.alpha = data.alpha_all_forward(data.fr);
                        data.dt = data.forward_dt;
                    end
                    
                    % Bottom function info
                    if m_iter == 1
                        noise_P = data.noise_guess*(rand-0.5);
                        B = B0 + data.p_guess(time)*(1+noise_P)*B1;
                        p_guess_vec(data.fr) = data.p_guess(time)*(1+noise_P);
                    else
                        B = B0 + p_guess(data.fr)*B1; 
                    end
                    B_quad = (vander_pd(:,1:data.pd)*B')';
                    
                    if data.fr == 1
                        B_all(:,:,1) = B(2:ex+1,:);
                    end

                    % RK step 1
                    [H,Q] = SWE_Update.runge_kutta(H0,Q0,0,0,B,B_quad,data,rk_step);
                    if data.accuracy_test == 0
                        [H,Q,~] = SWE_Limiters.slope_lim_characteristic(H,Q,B,data);
                    end

                    B0 = B0(2:end-1,:);
                    B1 = B1(2:end-1,:);
                else
            
                    % Bottom function info
                    if m_iter == 1
                        noise_P = data.noise_guess*(rand-0.5);
                        if rk_step == 2
                            B = B0 + data.p_guess(time+data.dt)*(1+noise_P)*B1;
                        elseif rk_step == 3  
                            B = B0 + data.p_guess(time+0.5*data.dt)*(1+noise_P)*B1;
                        end
                    else
                        if rk_step == 2 
                            B = B0 + p_guess(data.fr+1)*B1;
                        elseif rk_step == 3
                            % Use quadratic interpolation to get p(t+0.5*dt)
                            if data.fr > 1
                                p_rk3 = (-1/8)*p_guess(data.fr-1) + (3/4)*p_guess(data.fr) + (3/8)*p_guess(data.fr+1);
                            else
                                p_rk3 = (3/8)*p_guess(data.fr) + (3/4)*p_guess(data.fr+1) - (1/8)*p_guess(data.fr+2);
                            end
                            B = B0 + p_rk3*B1; 
                        end
                    end
                    B_quad = (vander_pd(:,1:data.pd)*B')';
                    
                    % B.C. padding
                    H0 = SWE_SetUp.bc_padding(H,bc_type);
                    Q0 = SWE_SetUp.bc_padding(Q,bc_type);
                    B  = SWE_SetUp.bc_padding(B,bc_type);

                    % RK steps 2 & 3
                    [H,Q] = SWE_Update.runge_kutta(H0,Q0,Hn,Qn,B,B_quad,data,rk_step);
                    if data.accuracy_test == 0
                        [H,Q,~] = SWE_Limiters.slope_lim_characteristic(H,Q,B,data);
                    end
                end

                % Clean Up B
                B = B(2:end-1,:);
            end


            %% ---------------------------------------------------------------------------
            %     Save Data as Needed
            % ----------------------------------------------------------------------------

            % Save variables for all time and space
            H_all(:,:,data.fr+1) = H;
            Q_all(:,:,data.fr+1) = Q;
            B_all(:,:,data.fr+1) = B;
            if m_iter == 1
                dt_all_forward(data.fr) = data.dt;
                alpha_all_forward(data.fr) = data.alpha;
            end

            
            %% ---------------------------------------------------------------------------
            %     Advance Time 
            % ----------------------------------------------------------------------------

            % Update time step counter
            data.fr = data.fr + 1;
            timestep = timestep + 1;
            
            % Update current time
            time = time + data.dt;
           

        end
    end

    if data.accuracy_test == 0
        % End timer & print time statement
        program_time = toc;
        if mod(m_iter,10) == 0
            disp(['Iteration ' num2str(m_iter) ' run time with ', num2str(ex), ' cells: ', num2str(program_time)]);
            disp(' ');
        end
    end
    
    % Remove extra 0s if necessary
    H_all = H_all(:,:,1:data.fr);
    Q_all = Q_all(:,:,1:data.fr);
    B_all = B_all(:,:,1:data.fr);
    
    filename_str = data.filename_str;
        
    if (m_iter == 1) 
        if data.accuracy_test == 0
            % Convert noisy data to more coarse time scale
            t_steps_forward = cumsum(dt_all_forward(1:data.fr-1));
            t_steps_forward = [0, t_steps_forward];
            data.t_steps_forward = t_steps_forward;

            t_steps_meas = data.t_steps_meas;
            t_steps_meas_idx = zeros(length(t_steps_forward),1);
            for t = 1:length(t_steps_forward)
                [~,t_steps_meas_idx(t)] = min(abs(t_steps_forward(t)-t_steps_meas));
            end

            data.H_noisy = data.H_noisy(:,:,t_steps_meas_idx);
            data.Q_noisy = data.Q_noisy(:,:,t_steps_meas_idx);
        end

        % Save Info
        data.dt_all_forward = dt_all_forward(1:data.fr-1);
        data.alpha_all_forward = alpha_all_forward(1:data.fr-1);
        
        if data.accuracy_test == 0
            noise_P = data.noise_guess*(rand-0.5);
            p_guess_vec(data.fr) = data.p_guess(time)*(1+noise_P);
            data.p_guess = p_guess_vec;

            file_pguess = sprintf('%s/SWE_PGuessNoisy_%s.csv',data.front_path,filename_str);
            file_dt_forward = sprintf('%s/SWE_ForwardTimeSteps_%s.csv',data.front_path,filename_str);
            file_dt_meas_idx = sprintf('%s/SWE_MeasTimeAdjIdx_%s.csv',data.front_path,filename_str);

            dlmwrite(file_pguess, data.p_guess,'delimiter',',','precision',20)
            dlmwrite(file_dt_forward, t_steps_forward,'delimiter',',','precision',20)
            dlmwrite(file_dt_meas_idx, t_steps_meas_idx,'delimiter',',','precision',20)

            % Convert noisy data to more coarse spatial mesh
            [data.H_noisy,~] = SWE_SetUp.convert_coarse(data.H_noisy,data,'forward',2); 
            [data.Q_noisy,~] = SWE_SetUp.convert_coarse(data.Q_noisy,data,'forward',2); 
        end
    end

    if data.accuracy_test == 0
        % Calculate residue
        H_res_all = H_all - data.H_noisy(2:ex+1,:,:);
        Q_res_all = Q_all - data.Q_noisy(2:ex+1,:,:); 

        if size(H_res_all,2) == 1 || size(H_res_all,2) == 2
            H_res_ends(:,1) = H_res_all(1,1,:);
            H_res_ends(:,2) = H_res_all(end,1,:);
            Q_res_ends(:,1) = Q_res_all(1,1,:);
            Q_res_ends(:,2) = Q_res_all(end,1,:);
        elseif size(H_res_all,2) == 3 || size(H_res_all,2) == 4
            H_res_ends(:,1) = H_res_all(1,1,:) + (1/12)*H_res_all(1,3,:);
            H_res_ends(:,2) = H_res_all(end,1,:) + (1/12)*H_res_all(end,3,:);
            Q_res_ends(:,1) = Q_res_all(1,1,:) + (1/12)*Q_res_all(1,3,:);
            Q_res_ends(:,2) = Q_res_all(end,1,:)+ (1/12)*Q_res_all(end,3,:);
        end

        % Cost function
        data.cost_func = 0.5*data.dt*sum(sum((H_res_ends(:,1)).^2,2)) + ...
            0.5*data.dt*sum(sum((Q_res_ends(:,1)).^2,2)) + ...
            0.5*data.dt*sum(sum((H_res_ends(:,2)).^2,2)) + ...
            0.5*data.dt*sum(sum((Q_res_ends(:,2)).^2,2));

        % Adjoint problem
        data = SWE_Inverse.adjoint_dg(H_all,Q_all,B_all,H_res_all,Q_res_all,data);
    end

    
    %% -----------------------------------------------------------------------------------
    %     Save Data for Error Calculations
    % ------------------------------------------------------------------------------------
        
    % Save predicted data for all time and space
    fileH_pred = sprintf('%s/SWE_PredData_%s_H.csv',data.front_path,filename_str);
    fileQ_pred = sprintf('%s/SWE_PredData_%s_Q.csv',data.front_path,filename_str);
    fileB_pred = sprintf('%s/SWE_PredData_%s_B.csv',data.front_path,filename_str);

    dlmwrite(fileH_pred, H_all,'delimiter',',','precision',20)
    dlmwrite(fileQ_pred, Q_all,'delimiter',',','precision',20)
    dlmwrite(fileB_pred, B_all,'delimiter',',','precision',20)

    if data.accuracy_test == 0
        % Save residue for all time and space
        fileH_res = sprintf('%s/SWE_ResData_%s_H.csv',data.front_path,filename_str);
        fileQ_res = sprintf('%s/SWE_ResData_%s_Q.csv',data.front_path,filename_str);

        dlmwrite(fileH_res, H_res_all,'delimiter',',','precision',20)
        dlmwrite(fileQ_res, Q_res_all,'delimiter',',','precision',20)
    end
    
    % Save predicted data at final time
    fileH_T = sprintf('%s/SWE_PredData_%s_H_FinalTime.csv',data.front_path,filename_str);
    fileQ_T = sprintf('%s/SWE_PredData_%s_Q_FinalTime.csv',data.front_path,filename_str);
    fileB_T = sprintf('%s/SWE_PredData_%s_B_FinalTime.csv',data.front_path,filename_str);

    dlmwrite(fileH_T, H,'delimiter',',','precision',20)
    dlmwrite(fileQ_T, Q,'delimiter',',','precision',20)
    dlmwrite(fileB_T, B,'delimiter',',','precision',20)

    data.H_all_pred_current = H_all;
    data.Q_all_pred_current = Q_all;
    data.B_all_pred_current = B_all;
    
    if data.accuracy_test == 0
        % Save Adjoint solution data
        fileS1_T = sprintf('%s/SWE_AdjData_%s_S1_FinalTime.csv',data.front_path,filename_str);
        fileS2_T = sprintf('%s/SWE_AdjData_%s_S2_FinalTime.csv',data.front_path,filename_str);

        dlmwrite(fileS1_T, data.S1,'delimiter',',','precision',20)
        dlmwrite(fileS2_T, data.S2,'delimiter',',','precision',20)
    end
end
    