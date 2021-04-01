classdef SWE_Inverse
    methods(Static)

        %% -----------------------------------------------------------------------
        %     Update p^k
        % ------------------------------------------------------------------------
                
        function [data, Z, m_iter] = inverse(Z,m_iter,data)

            % Unpack variables
            ex = data.forward_ex;
            pd = data.pd;
            xquad = data.xquad;
            wquad = data.wquad;
            coarse_ex = data.adjoint_ex;
            B1 = data.B1;
            lr = data.lr;
            lr_type = data.lr_type;
            alpha = lr;
            gammaL = data.gammaL;
            gammaH = data.gammaH;
            lambda = data.lambda;
            P0 = data.p0;
            P = data.p_guess;
            del_Lambda_star = data.del_Lambda_star;

            % Convert B1 to coarse mesh
            ratio = ex/coarse_ex;
            if ratio == 1
                B1_coarse = B1;
            else
                num_quad = ceil(pd/ratio);
                xquad_c = lgwt(num_quad,-0.5,0.5);
                if mod(ratio,2) == 0
                    shift = 0.5*(1/ratio)*(-(ratio-1):2:ratio-1);
                else
                    shift = (1/ratio)*(-0.5*(ratio-1):0.5*(ratio-1));
                end
                x_coarse = reshape(xquad_c/ratio + shift,[],1);
                vander_adapt = xquad_c.^[0:pd-1];
                B1_quad = reshape(vander_adapt*B1',[],coarse_ex)';
                vander_coarse = x_coarse.^[0:pd-1];
                B1_coarse = (vander_coarse\B1_quad')';
            end

            % Calculate gradient of J_0 = Integral_a^b (b_1)_x \delta \Lambda*
            B1_prime = 0;
            for i = 1:pd-1
                B1_prime = B1_prime + i*B1_coarse(:,i+1).*((xquad').^(i-1));
            end
            grad_J0 = zeros(size(P));
            for et = 1:size(del_Lambda_star,1)
                grad_J0(et) = sum(sum(B1_prime'.*squeeze(del_Lambda_star(et,:,:))'.*wquad));

            end

            % Proximal gradient descent step size
            [alpha, data] = SWE_SetUp.lr_adjust(lr,lr_type,grad_J0,m_iter,data);

            % Use proximal gradient descent to update p^m
            omega = 2*P - Z - alpha.*grad_J0;
            
            % Periodic Laplacian
            % Finite difference 2nd derivative: (u^{n} - 2 u^{n-1} + u^{n-2})/dt^2
            Lap = diag( ones(1,length(P)-1) , 1 ) - 2*diag( ones(1,length(P)) ,0 ) + diag( ones(1,length(P)-1) , -1 ) ; 
            Lap(1,end) = 1;
            Lap(end,1) = 1;
            Lap = (data.dt).^(-2)*Lap; 
                                 
            omega = transpose( ( eye(length(P)) - data.g_hat*gammaH*Lap ) \ transpose(omega) ) ;
            Z = Z + lambda*(omega - P);
            pnew = sign(Z-P0).*max(abs(Z-P0)-alpha*data.g_hat*gammaL,0)+P0; % add back the 1

            % Save error
            data.iter_err1(m_iter) = data.cost_func;
            reg_L = sum(abs(pnew-P0));
            reg_H = sum((([pnew(2:end) pnew(1)] - [pnew(end) pnew(1:end-1)])./ (2*data.dt)).^2);  
            data.iter_err_reg(m_iter) = data.cost_func + (gammaL*reg_L + gammaH*reg_H);
            if data.iter_err1(m_iter) < data.best_err1
                % Save info on best iteration
                data.best_err1 = data.iter_err1(m_iter);
                data.best_iter1 = m_iter;
                data.best_p1 = pnew;
                data.H_all_pred_best1 = data.H_all_pred_current;
                data.Q_all_pred_best1 = data.Q_all_pred_current;
                data.B_all_pred_best1 = data.B_all_pred_current;
                
                % Also calculate the regularization values
                reg_L = sum(abs(pnew-P0));
                reg_H = sum((([pnew(2:end) pnew(1)] - [pnew(end) pnew(1:end-1)])./ (2*data.dt)).^2);
                data.reg = (gammaL*reg_L + gammaH*reg_H);
                
            end
            data.best_err_all1(m_iter) = data.best_err1;

            if data.plotting_on == 1
                % Plot & save results every 'plot_times' iterations
                plot_times = data.plot_iters;
                if mod(m_iter,plot_times) == 0 || m_iter == 1
                    SWE_PostProcess.plot_pIteratations(m_iter,plot_times,P,pnew,data) 
                    if m_iter > 1
                        disp('')
                    end
                end
            end
            
            data.p_all(m_iter,:) = P;
            file_p_all = sprintf('%s/SWE_PData_%s.csv',data.front_path,data.filename_str);
            dlmwrite(file_p_all, data.p_all,'delimiter',',','precision',20)
                
            % Update guess
            P = pnew;
            data.p_guess = P;
            
            % Update iteration Counter
            m_iter = m_iter + 1;

        end
        
        
        %% ----------------------------------------------------------------
        %     DG Driver for Adjoint Problem
        % -----------------------------------------------------------------
                
        function data = adjoint_dg(H,Q,B,H_res,Q_res,data)
            % Unpack Variables
            dt_all = data.dt_all_forward; % assume same timestep size for coarser mesh
            alpha_all = data.alpha_all_forward;
            g = data.g;
            thresh = data.thresh;
            redef_method = data.redef_method;
            coarse_ex = data.adjoint_ex;
            coarse_pd = data.adjoint_k+1;
            vander_pd = data.vander_pd;
            
            % Coarse Mesh Set Up
            [H_c, H_CA_c] = SWE_SetUp.convert_coarse(H, data, 'adjoint', 2);
            [Q_c, Q_CA_c] = SWE_SetUp.convert_coarse(Q, data, 'adjoint', 2);
            [B_c, ~] = SWE_SetUp.convert_coarse(B, data, 'adjoint', 2);
            [~, H_res_CA_c] = SWE_SetUp.convert_coarse(H_res, data, 'adjoint', 1);
            [~, Q_res_CA_c] = SWE_SetUp.convert_coarse(Q_res, data, 'adjoint', 1);
            U_CA_c = zeros(size(H_CA_c));
            for t = 1:length(dt_all)
                [~,U_CA_c(t,:)] = SWE_SetUp.vel_redefine(squeeze(H_CA_c(t,:,:)),squeeze(Q_CA_c(t,:,:)),thresh,redef_method);
            end
            
            % Boundary Conditions for S1 & S2
            Sig_L = zeros(length(dt_all)+1,2);
            Sig_R = zeros(length(dt_all)+1,2);
            for t = length(dt_all)+1:-1:1
                A_transpose_L = [0, g*H_CA_c(t,1)-U_CA_c(t,1)^2; 1, 2*U_CA_c(t,1)];
                A_transpose_R = [0, g*H_CA_c(t,end)-U_CA_c(t,end)^2; 1, 2*U_CA_c(t,end)];
                Sig_L(end-t+1,:) = A_transpose_L\[-H_res_CA_c(t,1); -Q_res_CA_c(t,1)];
                Sig_R(end-t+1,:) = A_transpose_R\[H_res_CA_c(t,end); Q_res_CA_c(t,end)];
            end
            
            % Initial Conditions
            S1 = zeros(coarse_ex,coarse_pd);
            S2 = zeros(coarse_ex,coarse_pd);
                    
            % Update in time going backwards
            for tstep = 1:length(dt_all)+1
                t_back = length(dt_all)+1 - tstep + 1;
                
                for rk_step = 1:3
                    
                    % Store polynomial on 1st rk steps
                    if rk_step == 1
                        S1n = S1;
                        S2n = S2;
                    end
                    
                    % Boundary Conditions
                    S1_L = [Sig_L(tstep,1), zeros(coarse_pd-1,1)];
                    S1_R = [Sig_R(tstep,1), zeros(coarse_pd-1,1)];
                    S2_L = [Sig_L(tstep,2), zeros(coarse_pd-1,1)];
                    S2_R = [Sig_R(tstep,2), zeros(coarse_pd-1,1)];
%                     S1_L = [Sig_L(tstep,1), zeros(1,coarse_pd-1)];
%                     S1_R = [Sig_R(tstep,1), zeros(1,coarse_pd-1)];
%                     S2_L = [Sig_L(tstep,2), zeros(1,coarse_pd-1)];
%                     S2_R = [Sig_R(tstep,2), zeros(1,coarse_pd-1)];
                    
                    S1 = [S1_L; S1; S1_R];
                    S2 = [S2_L; S2; S2_R];
                    
                    % Time step and LF value
                    dt = dt_all(end);
                    alpha = alpha_all(end);
                    
                    % RK update
                    [S1,S2] = SWE_Update.runge_kutta_adjoint(S1,S2,S1n,S2n,squeeze(H_c(:,:,t_back)),squeeze(Q_c(:,:,t_back)),squeeze(B_c(:,:,t_back)),dt,alpha,data,rk_step);
                end
                
                data.S1 = S1;
                data.S2 = S2;
                
                S1_all(t_back,:,:) = S1;
                S2_all(t_back,:,:) = S2;
            end
            
            % For the gradient calculation
            S2_all_quad = zeros(size(S2_all,1),size(S2_all,2),size(vander_pd,1));
            H_c_quad = zeros(size(S2_all,1),size(S2_all,2),size(vander_pd,1));
            for t = 1:size(S2_all,1)
                % Evaluate s2 & H_c at quad points
                S2_all_quad(t,:,:) = (vander_pd(:,1:coarse_pd)*squeeze(S2_all(t,:,:))')';
                H_c_quad(t,:,:) = (vander_pd(:,1:coarse_pd)*squeeze(H_c(2:coarse_ex+1,:,t))')';    
            end
            data.del_Lambda_star = -S2_all_quad.*data.g.*H_c_quad;

        end


    end
end
