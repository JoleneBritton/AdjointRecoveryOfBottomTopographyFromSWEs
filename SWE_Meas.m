function data = SWE_Meas(data)

    %% -----------------------------------------------------------------------------------
    %     Parameters & Initial Conditions
    % ------------------------------------------------------------------------------------
    
    % Unpack variables
    ex = data.meas_ex;
    dx = data.meas_dx;
    x_center = data.meas_mesh_center;
    pd = data.meas_k+1;

    data.pd = pd;
    data.ex = ex;
    data.dx = dx;
    data.x_center = x_center;
    bc_type = data.bc_type;
    vander_pd = data.vander_pd;
    
    % Information for P
    p_true = data.p_true;

    % Initial conditions for H,Q & Determine B
    [H,Q,B0,B1] = SWE_SetUp.initial_conditions(data,pd,ex,dx,x_center);
    data.B0 = B0;
    data.B1 = B1;

    et = 500;

    % Initialize matrices for holding all data for space and time
    H_all = zeros(ex,data.pd,et);
    [Q_all, B_all] = deal(H_all);
    dt_all_meas = zeros(et,1);
    alpha_all_meas = zeros(et,1);
    
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
        
        for timestep = 1:data.meas_nt

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

                %% -----------------------------------------------------------------------
                %     Find dt & Runge-Kutta Steps
                % ------------------------------------------------------------------------
                    
                if rk_step == 1

                    % Global Lax-Friedrichs coefficient
                    data.alpha = SWE_SetUp.swe_alpha(H0,Q0,data);

                    % Determine dt
                    data.dt = data.meas_dt;
                    
                    % Bottom function info
                    B = B0 + p_true(time)*B1;
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
                    if rk_step == 2
                        B = B0 + p_true(time + data.dt)*B1;
                    elseif rk_step == 3
                        B = B0 + p_true(time + 0.5*data.dt)*B1;
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
            dt_all_meas(data.fr) = data.dt;
            alpha_all_meas(data.fr) = data.alpha;

            
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

    % End timer
    program_time = toc;
    disp(['Measured data computation run time with ', num2str(ex), ' cells: ', num2str(program_time)]);
    disp(' ');

    % Remove extra 0s if necessary
    H_all = H_all(:,:,1:data.fr);
    Q_all = Q_all(:,:,1:data.fr);
    B_all = B_all(:,:,1:data.fr);
    dt_all_meas = dt_all_meas(1:data.fr-1);
    alpha_all_meas = alpha_all_meas(1:data.fr-1);
        
    % Make noisy measured data
    noise_H = data.noise_meas*(rand(size(H_all))-0.5);
    noise_Q = data.noise_meas*(rand(size(Q_all))-0.5);
    data.H_noisy = H_all.*(1+noise_H);
    data.Q_noisy = Q_all.*(1+noise_Q);
    
    
    %% -----------------------------------------------------------------------------------
    %     Save Data for Error Calculations
    % ------------------------------------------------------------------------------------

    % Write all forward solutions to class instance
    data.H_all_meas = H_all;
    data.Q_all_meas = Q_all;
    data.B_all_meas = B_all;

    file_H = sprintf('%s/SWE_Meas_%s_H.csv',data.front_path,data.filename_str);
    file_Q = sprintf('%s/SWE_Meas_%s_Q.csv',data.front_path,data.filename_str);
    file_B = sprintf('%s/SWE_Meas_%s_B.csv',data.front_path,data.filename_str);

    dlmwrite(file_H, H_all,'delimiter',',','precision',20)
    dlmwrite(file_Q, Q_all,'delimiter',',','precision',20)
    dlmwrite(file_B, B_all,'delimiter',',','precision',20)

    % Save info about time steps
    t_steps_meas = cumsum(dt_all_meas);
    t_steps_meas = [0; t_steps_meas];
    data.t_steps_meas = t_steps_meas;
    data.dt_all_meas = dt_all_meas;
    data.alpha_all_meas = alpha_all_meas;
    
    file_tsteps = sprintf('%s/SWE_tsteps_meas_%s.csv',data.front_path,data.filename_str);
    file_dt_all = sprintf('%s/SWE_dt_meas_%s.csv',data.front_path,data.filename_str);
    file_alpha_all = sprintf('%s/SWE_alpha_meas_%s.csv',data.front_path,data.filename_str);

    dlmwrite(file_tsteps, t_steps_meas,'delimiter',',','precision',20)
    dlmwrite(file_dt_all, dt_all_meas,'delimiter',',','precision',20)
    dlmwrite(file_alpha_all, alpha_all_meas,'delimiter',',','precision',20)

end
    
