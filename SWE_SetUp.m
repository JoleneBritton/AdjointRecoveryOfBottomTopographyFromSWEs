
classdef SWE_SetUp
    methods(Static)
            
        % --------------------------------------------------------------------------------
        % P^k IC Projection function
        % --------------------------------------------------------------------------------
        
        function [H,Q,B0,B1] = initial_conditions(data,pd,ex,dx,x_center)
            
            % Unpack variables
            xquad = data.xquad;
            wquad = data.wquad;
                     
            % Mass Matrix  
            M = dx*SWE_SetUp.mass_matrix(pd);

            % Initialize arrays
            icH = zeros(pd,ex);                % Initial condition for H
            [icQ, icB0, icB1] = deal(icH);     % Initial conditions for Q,Z,B

            Hf = zeros(pd,ex,length(xquad));    
            [Qf, B0f, B1f]= deal(Hf);
            
            % Get variable initial conditions
            Hfunc = data.h0_func;
            Qfunc = data.q0_func;
            B0func = data.b0_func;
            B1func = data.b1_func;
            
            % Polynomial Computation
            % Multiply by test function [1; x_i; x_i^2; ...]
            for i = 1:pd
                Hf(i,:,:) = (Hfunc(xquad,dx,x_center).*(xquad.^(i-1)))';
                Qf(i,:,:) = (Qfunc(xquad,dx,x_center).*(xquad.^(i-1)))';
                B0f(i,:,:) = (B0func(xquad,dx,x_center).*(xquad.^(i-1)))';
                B1f(i,:,:) = (B1func(xquad,dx,x_center).*(xquad.^(i-1)))';
            end

            % Intial condition function coefficients
            for i = 1:pd
                for j = 1:ex
                    icH(i,j) = dx*sum(squeeze(Hf(i,j,:)).*wquad);
                    icQ(i,j) = dx*sum(squeeze(Qf(i,j,:)).*wquad);
                    icB0(i,j) = dx*sum(squeeze(B0f(i,j,:)).*wquad);
                    icB1(i,j) = dx*sum(squeeze(B1f(i,j,:)).*wquad);
                end
            end

            % Initialize
            H = zeros(ex,pd);
            Q = zeros(ex,pd);
            B0 = zeros(ex,pd);
            B1 = zeros(ex,pd);

            % Determine polynomial coefficients
            for j = 1:ex
                H(j,:) = M\icH(:,j);
                Q(j,:) = M\icQ(:,j);
                B0(j,:) = M\icB0(:,j);
                B1(j,:) = M\icB1(:,j);
            end


        end
        
        % --------------------------------------------------------------------------------
        % Cell interface function
        % --------------------------------------------------------------------------------
        
        function U_CI = cell_interfaces_old(U,j,pd)
            
            U_CI(:,1) = U(:,j-1)'*(0.5.^[0:pd-1])';
            U_CI(:,2) = U(:,j)'*((-0.5).^[0:pd-1])';
            U_CI(:,3) = U(:,j)'*(0.5.^[0:pd-1])';
            U_CI(:,4) = U(:,j+1)'*((-0.5).^[0:pd-1])';
        end
        
        function U_CI = cell_interfaces(U,pd)
            if pd == 1
                U_CI(:,1) = U(1:end-1);
                U_CI(:,2) = U(2:end);
            else
                U_CI(:,1) = U(1:end-1,:)*(0.5.^[0:pd-1])';
                U_CI(:,2) = U(2:end,:)*((-0.5).^[0:pd-1])';
            end
        end
        
        
        % --------------------------------------------------------------------------------
        % Boundary conditions function
        % --------------------------------------------------------------------------------
        
        function U = bc_padding(U,bc_type)
            
            % Periodic B.C.
            if isequal(bc_type, 'periodic')
                U = [U(end,:); U; U(1,:)];

            % No B.C.
            elseif isequal(bc_type, 'none')
                U = [U(1,:); U; U(end,:)]; 

            end            
        end
        
        
        % --------------------------------------------------------------------------------
        % Mass Matrix
        % --------------------------------------------------------------------------------
        
        function M = mass_matrix(pd)

            M = zeros(pd);
            vec = zeros(2*pd-1,0);

            for i = 1:2:2*pd-1
                vec(i) = 1/(i*2^(i-1));
            end

            for i = 1:pd
                M(i,:) = vec(i:pd+i-1);
            end 
        end
        
        
        % --------------------------------------------------------------------------------
        % LF Coefficient Alpha function
        % --------------------------------------------------------------------------------
        
        function alpha = swe_alpha(H,Q,data)
            
            % Unpack variables
            g = data.g;
            pd = data.pd;
            ex = data.ex;
            thresh = data.thresh;
            redef_method = data.redef_method;
            
            % Values in transformed cell
            xx = linspace(-.5,.5,100);
            
            % Evaluate function at many points in each cell
            vander_pd = data.vander_pd;
            vander_pd = vander_pd(:,1:pd);            
            Happrox = vander_pd * H';
            Qapprox = vander_pd * Q';

            % Clean Up H,Q
            Happrox = Happrox(:,2:ex+1);
            Qapprox = Qapprox(:,2:ex+1);

            [~,Uapprox] = SWE_SetUp.vel_redefine(Happrox,Qapprox,thresh,redef_method);
            
            % Eigenvalue of Jacobian calcuation
            fprime1 = abs(Uapprox) + sqrt(g*max(0,Happrox));
            
            % Alpha is a global max
            alpha = 1.2 * max(max(fprime1));
        end
        
        
        % --------------------------------------------------------------------------------
        % Redefine velocity function
        % --------------------------------------------------------------------------------
        
        function [Q,U] = vel_redefine(H,Q,thresh,redef_method)
            if redef_method == 0
                U = Q./H;
                U(isnan(U)) = 0;
                U(H <= thresh) = 0;

            elseif redef_method == 1
                new_u = find((H - thresh^(1/4) < 0) & (H>0));
                if ~isempty(new_u)
                    disp('')
                end
                U = Q./H;   
                U(new_u) = U(new_u).*(10.^((2.*((sqrt(2)./sqrt(1+(thresh./H(new_u).^4)))-1)).^3));
                U(abs(H)<1e-12) = 0;
                U(isnan(U)) = 0;
            end

            Q = H.*U;
        end
        
        
        % --------------------------------------------------------------------------------
        % LR Scheduluer
        % --------------------------------------------------------------------------------
        
        function [alpha, data] = lr_adjust(lr,lr_type,grad_J0,m_iter,data)
            changes = 10;
            nums_time = floor(length(grad_J0)/changes);
            nums_iters = floor(data.num_of_iters/changes);
            
            if lr_type == 0 || nums_time == 0 || nums_iters == 0
                alpha = lr*ones(size(grad_J0));
                data.alpha_type = 'constant';
                if nums_time == 0 || nums_iters == 0
                    data.lr_type = 0;
                end
            elseif lr_type == 1
                dec_perct = 0.1;
                dec_amount = lr-lr*(1-dec_perct);
                for c = 1:changes+1
                    alpha(1+(c-1)*nums_time:c*nums_time) = lr.*(1-(c-1)*dec_perct);
                end
                alpha = alpha(1:length(grad_J0));
                data.alpha_type = ['decrease by ' num2str(dec_amount) ' every ' num2str(nums_time) ' timesteps'];
            elseif lr_type == 2
                dec_perct = 0.05;
                dec_amount = lr-lr*(1-dec_perct);
                mod_iter_index = ceil(m_iter*changes/data.num_of_iters);
                alpha = lr*(1 - (mod_iter_index-1) * dec_perct);
                data.alpha_type = ['decrease by ' num2str(dec_amount) ' every ' num2str(nums_iters) ' iterations'];
            elseif lr_type == 3
                inc_amount = 5;
                inc_freq = 10;
                max_lr = 20;
                if mod(m_iter,inc_freq) == 0
                    data.lr = min(data.lr + inc_amount, max_lr);
                end
                alpha = data.lr;
                data.alpha_type = ['increase by ' num2str(inc_amount) ' every ' num2str(inc_freq) ' iterations'];
            end
            data.lr_all(m_iter) = alpha(1);
        end
        
        
        % --------------------------------------------------------------------------------
        % Convert to Coarse Mesh
        % --------------------------------------------------------------------------------
        
        function [U_c, U_c_CA] = convert_coarse(U,data,coarse_type,j_start)
            
            % Extract variables
            if isequal(coarse_type, 'adjoint')
                dt_all = data.dt_all_forward; % assume adjoint uses same time mesh as forward
                coarse_pd = data.adjoint_k+1;
                fine_pd = data.forward_k+1;
                coarse_ex = data.adjoint_ex;
                fine_ex = data.forward_ex;
            elseif isequal(coarse_type, 'forward')
                dt_all = data.dt_all_forward;
                coarse_pd = data.forward_k+1;
                fine_pd = data.meas_k+1;
                coarse_ex = data.forward_ex;
                fine_ex = data.meas_ex;
            end
            bc_type = data.bc_type;
            xquad = data.xquad;
            wquad = data.wquad;
            vander_pd = fliplr(vander(xquad));
            
            ratio = fine_ex/coarse_ex; 
            
            %Initialize
            U_c = zeros(coarse_ex+2,coarse_pd,length(dt_all));
            U_c_quad = zeros(coarse_ex+2,length(xquad),length(dt_all));
            U_c_CA = zeros(length(dt_all),coarse_ex+2);
            
            if ratio > 1
                num_quad = ceil(coarse_pd/ratio);
                
                % Points to eval. fine function at to quadrature vals.
                xquad_c = lgwt(num_quad,-0.5,0.5);
                
                % Points used to get coarse coefficents
                if mod(ratio,2) == 0
                    shift = 0.5*(1/ratio)*(-(ratio-1):2:ratio-1);
                else
                    shift = (1/ratio)*(-0.5*(ratio-1):0.5*(ratio-1));
                end
                x_coarse = reshape(xquad_c/ratio + shift,[],1);
                
                % Use: To eval fine function at quadrature points
                % Size = (len(xquad_c)) x poly. deg. of fine function
                vander_adapt = xquad_c.^(0:fine_pd-1);
                
                % Use: In conversion to get coefficients for coarse function, 
                % Size = (len(x_coarse)) x poly. deg. of coarse function
                vander_coarse = x_coarse.^(0:coarse_pd-1);

                % Convert to coarse mesh
                for t = 1:length(dt_all)+1
                    U_quad = reshape(vander_adapt*squeeze(U(:,:,t))',[],coarse_ex);
                    U_c(:,:,t) = SWE_SetUp.bc_padding((vander_coarse\U_quad)',bc_type);
                    if coarse_pd == 1
                        U_c_quad(:,:,t) = vander_pd(:,1:coarse_pd)*squeeze(U_c(:,:,t))';
                    else
                        U_c_quad(:,:,t) = squeeze(U_c(:,:,t))*vander_pd(:,1:coarse_pd)';
                    end
                    U_c_CA(t,:) = squeeze(U_c_quad(:,:,t))*wquad;
                end
            elseif ratio == 1
                
                xquad_c = lgwt(coarse_pd, -0.5, 0.5);
                vander_fine = xquad_c.^(0:fine_pd-1);
                vander_coarse = xquad_c.^(0:coarse_pd-1);
                
                for t = 1:length(dt_all)+1
                    U_quad = reshape(vander_fine*squeeze(U(:,:,t))',[],coarse_ex);
                    U_c(:,:,t) = SWE_SetUp.bc_padding((vander_coarse\U_quad)',bc_type);
                    if isequal(coarse_type, 'adjoint')
                        U_c_quad(:,:,t) = (vander_pd(:,1:coarse_pd)*squeeze(U_c(:,:,t))')';
                    elseif isequal(coarse_type, 'forward')
                        U_c_quad(:,:,t) = (vander_pd(:,1:data.forward_k+1)*(squeeze(U_c(:,:,t))'))';
                    end
                    U_c_CA(t,:) = squeeze(U_c_quad(:,:,t))*wquad;
                end
            end
            disp('')
        end
    end
end
