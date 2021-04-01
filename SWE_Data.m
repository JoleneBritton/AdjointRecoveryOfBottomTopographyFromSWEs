classdef SWE_Data
    properties
        % Constants
        g;
        choice;
        thresh;
        redef_method;
        CFL;
        num_of_iters;
        plot_iters;
        plotting_on;
        
        % Hyperparameters
        lr;
        lr_type;
        lr_all;
        alpha_type;
        g_hat;
        gammaL;
        gammaH;
        lambda;
       
        % Domain Info
        pd;
        xquad;
        wquad;
        vander_pd;
        start_pt;
        end_pt;
        ex;
        dx;
        x_center;
        meas_mesh_center;
        meas_ex;
        meas_dx;
        meas_k;
        forward_mesh_center;
        forward_ex;
        forward_dx;
        forward_k;
        adjoint_mesh_center;
        adjoint_ex;
        adjoint_dx;
        adjoint_k
        B0;
        B1;
        
        % Time Info
        dt;
        et;
        fr;
        rk_step;
        T;
        t_steps_meas;
        t_steps_forward;
        dt_all_meas;
        dt_all_forward;
        alpha_all_meas;
        alpha_all_forward;
        meas_nt;
        meas_dt;
        forward_nt;
        forward_dt;
        adjoint_nt;
        adjoint_dt;
        
        % ICs
        h0_func;
        q0_func;
        b0_func;
        b1_func;
        
        % P Info
        p_true;
        p_guess;
        p0;
        noise_meas;
        noise_guess;
        
        % Problem classifiers
        bc_type;
        meas_data;
        
        % Plotting Info
        x_lim;
        y_lim;
        leg_loc;
        title_str;
        
        % Measured Data 
        H_all_meas;
        Q_all_meas;
        B_all_meas;
        HQbound_all_true;
        H_noisy;
        Q_noisy;
        H_noisy_bdry;
        Q_noisy_bdry;
        
        % Iteration Error
        cost_func;
        iter_err1;
        best_err1;
        best_iter1;
        best_p1;
        best_err_all1;
        reg;
        iter_err_reg;
        
        H_all_pred_current;
        Q_all_pred_current;
        B_all_pred_current;
        H_all_pred_best1;
        Q_all_pred_best1;
        B_all_pred_best1;
        p_all;
        
        % Other
        alpha;
        HQbound_all_pred;
        del_Lambda_star;
        filename_str;
        front_path;
        S1;
        S2;
        accuracy_test;

    end
    
    methods(Static)
        
        function data = data_params(k,nt,func_vec,meshes,hyparams,choice,accuracy_test)
            
            data = SWE_Data();
            
            % Constants
            data.g = 9.812;                      % Gravitational constant
            data.CFL = 1/10;                     % Multiplier for dt
            data.thresh = 1e-6;                  % Threshold value for redefining Q
            data.redef_method = 0;               % Method for redefining Q for small H
            data.num_of_iters = hyparams(6);              % Number of iterations for the forward/inverse problem
            data.plot_iters = hyparams(7);
            data.plotting_on = 1;
            data.accuracy_test = accuracy_test;
            
            % Inverse Problem Hyperparamaters
            data.lr = hyparams(1);                    % Proximal gradient learning rate
            data.lr_type = hyparams(2);               % Proximal gradient learning rate
            data.gammaL = hyparams(3);                % L1 Regularization constant
            data.gammaH = hyparams(4);                % H1 Regularization constant
            data.lambda = 1;                          % Relaxation parameter
            data.g_hat = hyparams(10);
            
            % Temporal Info
            data.T = hyparams(5);  
            
            % Spatial Info
            if func_vec(1) == 1
                % Oscillatory ICs, b_0 = 0
                data.start_pt = 0;
                data.end_pt = 2;
                data.h0_func = @(x,dx,x_center) (7 + exp(sin(20*pi*(x_center(1:end-1)+dx*x))));
                data.q0_func = @(x,dx,x_center) (sin(cos(20*pi*(x_center(1:end-1)+dx*x))));
                data.b1_func = @(x,dx,x_center) ((sin(0.5*pi*(x_center(1:end-1)+dx*x))).^2);
                data.b0_func = @(x,dx,x_center) (zeros(1,length(x_center)-1) + 0*x);
                data.bc_type = 'periodic';
            elseif func_vec(1) == 2
                % Trig ICs, b_0 != 0
                data.start_pt = 0;
                data.end_pt = 1;
                data.h0_func = @(x,dx,x_center) (7 + exp(sin(2*pi*(x_center(1:end-1)+dx*x))));
                data.q0_func = @(x,dx,x_center) cos(2*pi*(x_center(1:end-1)+dx*x));
                data.b1_func = @(x,dx,x_center) ((sin(pi*(x_center(1:end-1)+dx*x))).^2);
                data.b0_func = @(x,dx,x_center) cos(sin(2*pi*(x_center(1:end-1)+dx*x)));
                data.bc_type = 'periodic';
            end
            
            data.meas_k = k(1);
            data.forward_k = k(2);
            data.adjoint_k = k(3);
            
            data.meas_nt = nt(1);
            data.forward_nt = nt(2);
            data.adjoint_nt = nt(3);
            
            data.meas_dt = data.T / data.meas_nt;
            data.forward_dt = data.T / data.forward_nt;
            data.adjoint_dt = data.T / data.adjoint_nt;
            
            data.meas_ex = meshes(1); 
            data.forward_ex = meshes(2);              % Forward problem number of cells
            data.adjoint_ex = meshes(3);          % Adjoint problem number of cells
            
            data.meas_dx = (data.end_pt - data.start_pt)/data.meas_ex;
            data.forward_dx = (data.end_pt - data.start_pt)/data.forward_ex;
            data.adjoint_dx = (data.end_pt - data.start_pt)/data.adjoint_ex;
            
            data.meas_mesh_center = linspace(data.start_pt,data.end_pt,data.meas_ex+1) + 0.5*data.meas_dx;
            data.forward_mesh_center = linspace(data.start_pt,data.end_pt,data.forward_ex+1) + 0.5*data.forward_dx;
            data.adjoint_mesh_center = linspace(data.start_pt,data.end_pt,data.adjoint_ex+1) + 0.5*data.adjoint_dx;
            
            [data.xquad,data.wquad] = lgwt(5,-0.5,0.5); % Quadrature points and weights
            data.vander_pd = fliplr(vander(data.xquad));

            % P Functions & Noise
            data.p0 = 1;
            if data.T == 0.01
                p_true_mult = -50000;
            elseif data.T == 0.05
                p_true_mult = -10000;
            elseif data.T == 0.1
                p_true_mult = -1000;
            elseif data.T == 0.2
                p_true_mult = -700;
            elseif data.T == 0.3
                p_true_mult = -500;
            end
            %%% ----- True P values -----
            if func_vec(2) == 0
                % Constant Function
                p_true = @(t) 0*t + 1;
            elseif func_vec(2) == 1
                % 1 Bump - Centered
                p_true =  @(t) (exp(p_true_mult*(t-0.5*data.T).^2)+data.p0);
            elseif func_vec(2) == 2
                % 1 Bump - Left Side
                p_true =  @(t) (exp(p_true_mult*(t-(1/3)*data.T).^2)+data.p0);
            elseif func_vec(2) == 3 
                % 1 Bump - Right Side
                p_true =  @(t) (exp(p_true_mult*(t-(2/3)*data.T).^2)+data.p0);
            elseif func_vec(2) == 4
                % 2 Bumps - Same height and width
                p_true =  @(t) (exp(2*p_true_mult*(t-(1/4)*data.T).^2)+exp(2*p_true_mult*(t-(3/4)*data.T).^2)+data.p0);
            elseif func_vec(2) == 5
                % 2 Bumps - Right taller and skinnier
                p_true =  @(t) (exp(p_true_mult*(t-0.3*data.T).^2)+1.5*exp(2*p_true_mult*(t-0.7*data.T).^2)+data.p0);
            elseif func_vec(2) == 6
                % 2 Bumps - Left taller and fatter
                p_true =  @(t) (1.5*exp(p_true_mult*(t-0.3*data.T).^2)+exp(2*p_true_mult*(t-0.7*data.T).^2)+data.p0);
            elseif func_vec(2) == 7
                % 2 Bumps & 1 Well
                p_true =  @(t) (exp(4*p_true_mult*(t-0.25*data.T).^2)+1.5*exp(4*p_true_mult*(t-0.5*data.T).^2)-0.5*exp(4*p_true_mult*(t-0.75*data.T).^2)+data.p0);
            elseif func_vec(2) == 8
                p_true =  @(t) (t.^2 + 0.5*exp(-sin(2.5*pi/data.T*(t-0.01))));
            elseif func_vec(2) == 9
                % No temporal component
                p_true =  @(t) (0*t);
            elseif func_vec(2) == 10
                p_true =  @(t) t;
            elseif func_vec(2) == 11
                p_true =  @(t) 1000000*(t.^3);
            end

            %%% ----- Initial Guess for P -----
            if func_vec(3) == 0
                % Constant Function
                p_guess = @(t) 0*t;
            elseif func_vec(3) == 1
                % Constant Function
                p_guess = @(t) 0*t + 1;
            elseif func_vec(3) == 2
                % Oscillatory Function - Medium Amplitude
                const = (10*pi)/data.T;
                p_guess = @(t) 1.5*cos(const*t).^2 + 0.75;
            elseif func_vec(3) == 3
                % Oscillatory Function - Larger Amplitude
                const = (10*pi)/data.T;
                p_guess = @(t) 3*cos(const*t).^2 + 0.75;
            elseif func_vec(3) == 4
                % 1 Bump - Centered, Larger Amplitude
                p_guess = @(t) (1.5*exp(0.5*p_true_mult*(t-(1/2)*data.T).^2)+data.p0);
            elseif func_vec(3) == 5
                % 1 Bump - Left Side
                p_guess = @(t) (exp(p_true_mult*(t-(1/3)*data.T).^2)+data.p0);
            elseif func_vec(3) == 6
                % 1 Bump - Right Side
                p_guess = @(t) (exp(p_true_mult*(t-(2/3)*data.T).^2)+data.p0);
            elseif func_vec(3) == 7
                % sin^2 - Upside down
                p_guess = @(t) -2*(sin(pi*t/data.T)).^2+2;
            elseif func_vec(3) == 8
                % sin^2 - Large Amplitude
                p_guess = @(t) 4*(sin(pi*t/data.T)).^2;
            elseif func_vec(3) == 9
                % 2 Bumps - Right taller and skinnier
                p_guess =  @(t) (exp(p_true_mult*(t-0.3*data.T).^2)+1.5*exp(2*p_true_mult*(t-0.7*data.T).^2)+data.p0);
            elseif func_vec(3) == 10
                p_guess =  @(t) t;
            elseif func_vec(3) == 11
                p_guess =  @(t) 1000000*(t.^3);
            end
            data.p_true = p_true;
            data.p_guess = p_guess;
            data.noise_meas = hyparams(8);
            data.noise_guess = hyparams(9);
            
            % Classifiers
            data.choice = choice;                % Numerical example identifier
            filename = sprintf('Eq=%d_LF_lr=%g_gL=%g_gH=%g_ghat=%g_exF=%d_exC=%d_T=%g_nM=%g_nG=%g_itrs=%d',...
                choice,data.lr,data.gammaL,data.gammaH,data.g_hat,data.forward_ex,data.adjoint_ex,data.T,data.noise_meas,data.noise_guess,data.num_of_iters);
            data.filename_str = strrep(filename, '.', ',');
            
            data.front_path = sprintf('SWE_Eq%d/k=%d,%d,%d/%s',choice,k(1),k(2),k(3),data.filename_str);
        end
        
        
    end
    
    
end
