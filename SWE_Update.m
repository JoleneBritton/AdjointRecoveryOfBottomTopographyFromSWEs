classdef SWE_Update
    methods(Static)
        
        % Forward Runge-Kutta function
        function [Hnew, Qnew] = runge_kutta(H,Q,H0,Q0,B,B_quad,data,RKstep)
            

            %% ---------------------------------------------------------------------------------------
            %     Name Varargin Inputs
            % ----------------------------------------------------------------------------------------
           
            % Unpack variables
            g = data.g;
            alpha = data.alpha;
            pd = data.pd;
            ex = data.ex;
            dx = data.dx;
            dt = data.dt;
            xquad = data.xquad;
            wquad = data.wquad;
            thresh = data.thresh;
            redef_method = data.redef_method;
            
            % Cell interface values
            H_CI = SWE_SetUp.cell_interfaces(H,pd);
            Q_CI = SWE_SetUp.cell_interfaces(Q,pd);
            B_CI = SWE_SetUp.cell_interfaces(B,pd);
            
            % Values at quadrature points
            vander_pd = fliplr(vander(xquad));
            H_quad = (vander_pd(:,1:pd)*H')';
            Q_quad = (vander_pd(:,1:pd)*Q')';
            
            % Mass matrix
            M = dx.*SWE_SetUp.mass_matrix(pd);  

            % Test functions
            pwrs = 0:pd-1;
            prmpwrs = [0 0:pd-2];
            phi_rside  = 0.5.^pwrs';           % Test function for the left
            phi_lside = (-0.5).^pwrs';         % Test function for the right
            phi = (xquad.^pwrs)';
            phiPrime = (repmat(pwrs,[length(xquad),1]).*(xquad.^(prmpwrs)))';

            % LAR well balanced functions at interfaces
            B_star = max(B_CI(:,1),B_CI(:,2));
            H_star = max(0,H_CI + B_CI - B_star);
            Q_star = Q_CI;

            % Integral terms
            [fH_quad, fQ_quad] = SWE_Update.fluxes(H_quad(2:ex+1,:),Q_quad(2:ex+1,:),g,thresh,redef_method);
            
            H_integral = zeros(ex,pd);
            Q_integral = zeros(ex,pd);
            for i = 2:pd
                H_integral(:,i) = sum(fH_quad.*wquad'.*phiPrime(i,:),2);
                Q_integral(:,i) = sum(fQ_quad.*wquad'.*phiPrime(i,:),2);
            end
            
            % Flux terms
            [fH_star, fQ_star] = SWE_Update.fluxes(H_star,Q_star,g,thresh,redef_method); 
            [fH_CI, fQ_CI] = SWE_Update.fluxes(H_CI,Q_CI,g,thresh,redef_method);
            
            H_LF_flux_lside = 0.5*(fH_star(1:ex,1) + fH_star(1:ex,2) - alpha*(H_star(1:ex,2) - H_star(1:ex,1)));
            H_LF_flux_rside = 0.5*(fH_star(2:ex+1,1) + fH_star(2:ex+1,2) - alpha*(H_star(2:ex+1,2) - H_star(2:ex+1,1)));
            
            H_flux_lside = (H_LF_flux_lside + fH_CI(1:ex,2) - fH_star(1:ex,2)).*phi_lside';
            H_flux_rside = (H_LF_flux_rside + fH_CI(2:ex+1,1) - fH_star(2:ex+1,1)).*phi_rside';
            
            Q_LF_flux_lside = 0.5*(fQ_star(1:ex,1) + fQ_star(1:ex,2) - alpha*(Q_star(1:ex,2) - Q_star(1:ex,1)));
            Q_LF_flux_rside = 0.5*(fQ_star(2:ex+1,1) + fQ_star(2:ex+1,2) - alpha*(Q_star(2:ex+1,2) - Q_star(2:ex+1,1)));
            
            Q_flux_lside = (Q_LF_flux_lside + fQ_CI(1:ex,2) - fQ_star(1:ex,2)).*phi_lside';
            Q_flux_rside = (Q_LF_flux_rside + fQ_CI(2:ex+1,1) - fQ_star(2:ex+1,1)).*phi_rside';
            
            % Source term
            if pd == 1
                Bprime = zeros(size(H_quad(2:ex+1,:)));
            else
                Bprime = 0;
                for i = 1:pd-1
                    Bprime = Bprime + i*B(2:ex+1,i+1).*((xquad').^(i-1));
                end
            end
            
            source_quad = -g*H_quad(2:ex+1,:).*Bprime;
            
            Q_source = zeros(ex,pd);
            for i = 1:pd
                Q_source(:,i) = sum(source_quad.*wquad'.*phi(i,:),2);
            end
            
            % Update variables
            H_update = H_integral - H_flux_rside + H_flux_lside;
            Q_update = Q_integral - Q_flux_rside + Q_flux_lside + Q_source;

            Hnew = zeros(ex,pd);
            Qnew = zeros(ex,pd);

            for j = 1:ex 
                if RKstep == 1
                    Hnew(j,:) = H(j+1,:) + (M\(dt*H_update(j,:)'))';
                    Qnew(j,:) = Q(j+1,:) + (M\(dt*Q_update(j,:)'))';
                elseif RKstep == 2
                    Hnew(j,:) = 0.75*H0(j+1,:) + 0.25*H(j+1,:) + 0.25*(M\(dt*H_update(j,:)'))';
                    Qnew(j,:) = 0.75*Q0(j+1,:) + 0.25*Q(j+1,:) + 0.25*(M\(dt*Q_update(j,:)'))';
                else
                    Hnew(j,:) = (1/3)*H0(j+1,:) + (2/3)*H(j+1,:) + (2/3)*(M\(dt*H_update(j,:)'))';
                    Qnew(j,:) = (1/3)*Q0(j+1,:) + (2/3)*Q(j+1,:) + (2/3)*(M\(dt*Q_update(j,:)'))';
                end
            end
            
        end
        
        
        %% -------------------------------------------------------------------------------
        %     Flux Function
        % --------------------------------------------------------------------------------
        function [fA, fQ] = fluxes(H,Q,g,thresh,redef_method)

            [~,U] = SWE_SetUp.vel_redefine(H,Q,thresh,redef_method);

            % Calculate Flux Values at interfaces
            fA = Q;
            fQ = Q.*U + 0.5*g*H.^2;
        end
        
        
        %% -------------------------------------------------------------------------------
        % Adjoint Runge-Kutta function
        % --------------------------------------------------------------------------------
        function [S1new, S2new] = runge_kutta_adjoint(sig1,sig2,sig1n,sig2n,H,Q,B,dt,alpha,data,RKstep)
         
            % Unpack variables
            g = data.g;
            ex = data.adjoint_ex;
            xquad = data.xquad;
            wquad = data.wquad;
            thresh = data.thresh;
            redef_method = data.redef_method;
            coarse_dx = data.adjoint_dx;
            coarse_pd = data.adjoint_k+1;
            
            if data.accuracy_test == 1
                H = SWE_SetUp.bc_padding(H, data.bc_type);
                Q = SWE_SetUp.bc_padding(Q, data.bc_type);
                B = SWE_SetUp.bc_padding(B, data.bc_type);
            end
            
            % Mass matrix
            M = coarse_dx.*SWE_SetUp.mass_matrix(coarse_pd); 
            
            % Cell interface values
            S1_CI = SWE_SetUp.cell_interfaces(sig1,coarse_pd);
            S2_CI = SWE_SetUp.cell_interfaces(sig2,coarse_pd);
            H_CI = SWE_SetUp.cell_interfaces(H,coarse_pd);
            Q_CI = SWE_SetUp.cell_interfaces(Q,coarse_pd);
            
            % Values at quadrature points
            vander_pd = fliplr(vander(xquad));
            S1_quad = (vander_pd(:,1:coarse_pd)*sig1')';
            S2_quad = (vander_pd(:,1:coarse_pd)*sig2')';
            H_quad = (vander_pd(:,1:coarse_pd)*H')';
            Q_quad = (vander_pd(:,1:coarse_pd)*Q')';
            [~,U_quad] = SWE_SetUp.vel_redefine(H_quad,Q_quad,thresh,redef_method);
            
            % Polynomial coefficients for U
            for i = 1:coarse_pd
                Uf(i,:,:) = (U_quad'.*(xquad.^(i-1)))';
            end
            for i = 1:coarse_pd 
                for j = 1:ex
                    icU(i,j) = coarse_dx*sum(squeeze(Uf(i,j,:)).*wquad);
                end
            end
            U = zeros(coarse_pd,ex);
            for j = 1:ex
                U(:,j) = M\icU(:,j);
            end
            U = SWE_SetUp.bc_padding(U',data.bc_type);

            % Test functions
            pwrs = 0:coarse_pd-1;
            prmpwrs = [0 0:coarse_pd-2];
            phi_rside  = 0.5.^pwrs';           % Test function for the left
            phi_lside = (-0.5).^pwrs';         % Test function for the right
            phi = (xquad.^pwrs)';
            phiPrime = (repmat(pwrs,[length(xquad),1]).*(xquad.^(prmpwrs)))';

            % Integral terms
            [fS1_quad, fS2_quad] = SWE_Update.adjoint_fluxes(S1_quad(2:ex+1,:),S2_quad(2:ex+1,:),H_quad(2:ex+1,:),Q_quad(2:ex+1,:),g,thresh,redef_method);
            
            S1_integral = zeros(ex,coarse_pd);
            S2_integral = zeros(ex,coarse_pd);
            if coarse_pd > 1
                for i = 2:coarse_pd
                    S1_integral(:,i) = sum(fS1_quad.*wquad'.*phiPrime(i,:),2);
                    S2_integral(:,i) = sum(fS2_quad.*wquad'.*phiPrime(i,:),2);
                end
            end
            
            % Flux terms
            [fS1_CI, fS2_CI] = SWE_Update.adjoint_fluxes(S1_CI,S2_CI,H_CI,Q_CI,g,thresh,redef_method);
             
            S1_flux_lside = 0.5*(fS1_CI(1:ex,1) + fS1_CI(1:ex,2) - alpha*(S1_CI(1:ex,2) - S1_CI(1:ex,1))).*phi_lside';
            S1_flux_rside = 0.5*(fS1_CI(2:ex+1,1) + fS1_CI(2:ex+1,2) - alpha*(S1_CI(2:ex+1,2) - S1_CI(2:ex+1,1))).*phi_rside';
             
            S2_flux_lside = 0.5*(fS2_CI(1:ex,1) + fS2_CI(1:ex,2) - alpha*(S2_CI(1:ex,2) - S2_CI(1:ex,1))).*phi_lside';
            S2_flux_rside = 0.5*(fS2_CI(2:ex+1,1) + fS2_CI(2:ex+1,2) - alpha*(S2_CI(2:ex+1,2) - S2_CI(2:ex+1,1))).*phi_rside';
             
            % Derivatives for source terms
            if coarse_pd == 1
                Hprime = 0.5*(H(3:end)-H(1:end-2)).*ones(ex,5);
                Uprime = 0.5*(U(3:end)-U(1:end-2)).*ones(ex,5);
                Bprime = 0.5*(B(3:end)-B(1:end-2)).*ones(ex,5);
            else
                Hprime = 0;
                Uprime = 0;
                Bprime = 0;
                for i = 1:coarse_pd-1
                    Hprime = Hprime + i*H(2:ex+1,i+1).*((xquad').^(i-1));
                    Uprime = Uprime + i*U(2:ex+1,i+1).*((xquad').^(i-1));
                    Bprime = Bprime + i*B(2:ex+1,i+1).*((xquad').^(i-1));
                end
            end
            
            % Source terms
            S1_source_quad = (g*(Hprime+Bprime) - 2*U_quad(2:ex+1,:).*Uprime).*S2_quad(2:ex+1,:);
            S2_source_quad = 2*Uprime.*S2_quad(2:ex+1,:);
            
            S1_source = zeros(ex,coarse_pd);
            S2_source = zeros(ex,coarse_pd);
            for i = 1:coarse_pd
                S1_source(:,i) = sum(S1_source_quad.*wquad'.*phi(i,:),2);
                S2_source(:,i) = sum(S2_source_quad.*wquad'.*phi(i,:),2);
            end
            
            % Update variables
            S1_update = S1_integral - S1_flux_rside + S1_flux_lside + S1_source;
            S2_update = S2_integral - S2_flux_rside + S2_flux_lside + S2_source;

            S1new = zeros(ex,coarse_pd);
            S2new = zeros(ex,coarse_pd);

            % RK3 calculations
            for j = 1:ex 
                if RKstep == 1
                    S1new(j,:) = sig1(j+1,:) + (M\(dt*S1_update(j,:)'))';
                    S2new(j,:) = sig2(j+1,:) + (M\(dt*S2_update(j,:)'))';
                elseif RKstep == 2
                    S1new(j,:) = 0.75*sig1n(j,:) + 0.25*sig1(j+1,:) + 0.25*(M\(dt*S1_update(j,:)'))';
                    S2new(j,:) = 0.75*sig2n(j,:) + 0.25*sig2(j+1,:) + 0.25*(M\(dt*S2_update(j,:)'))';
                else
                    S1new(j,:) = (1/3)*sig1n(j,:) + (2/3)*sig1(j+1,:) + (2/3)*(M\(dt*S1_update(j,:)'))';
                    S2new(j,:) = (1/3)*sig2n(j,:) + (2/3)*sig2(j+1,:) + (2/3)*(M\(dt*S2_update(j,:)'))';
                end
            end
            
        end
        
        
        %% -------------------------------------------------------------------------------
        %     Adjoint Flux Function
        % --------------------------------------------------------------------------------
        function [fs1, fs2] = adjoint_fluxes(S1,S2,H,Q,g,thresh,redef_method)

            [~,U] = SWE_SetUp.vel_redefine(H,Q,thresh,redef_method);

            % Calculate Flux Values at interfaces
            fs1 = (g*H-U.^2).*S2;
            fs2 = S1 + (2*U).*S2;
        end
        
        
    end
end
