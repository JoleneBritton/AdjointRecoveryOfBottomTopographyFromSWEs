classdef SWE_Limiters
    methods(Static)          
        
        %% -------------------------------------------------------------------------------
        %     Slope Limiter - Check if flagges along characteristic direction
        % --------------------------------------------------------------------------------
        
        function [H,Q,limiter_on] = slope_lim_characteristic(H,Q,B,data)
            
            % Unpack variables
            g = data.g;
            ex = data.ex;
            pd = data.pd;
            bc_type = data.bc_type;
            dx = data.dx;
            N = ex;
            m = pd-1;
            
            if data.meas_data == 1
%                 lim_type = 'weno';
                lim_type = 'minmod';
            else
%                 lim_type = 'minmod';
                lim_type = 'weno';
            end
            
            % Convert to nodal DG points
            if pd == 2
                x_gl = [-0.5;0.5];
            elseif pd == 3
                x_gl = [-0.5;0;0.5];
            elseif pd == 4
                x_gl = [-0.5;-sqrt(5)/10;sqrt(5)/10;0.5];
            end
            V_gl = fliplr(vander(x_gl));
            
            H = SWE_SetUp.bc_padding(H,bc_type)';
            Q = SWE_SetUp.bc_padding(Q,bc_type)';
            HB = H+B';
            
            VX = (data.end_pt - data.start_pt)*(0:N)/N + data.start_pt;
            x_gl_all = ones(m+1,1)*VX(1:N)+0.5*(2*x_gl+1)*(VX(2:N+1)-VX(1:N));
            
            HB_gl = (V_gl(:,1:pd)*HB(:,2:ex+1));
            H_gl = (V_gl(:,1:pd)*H(:,2:ex+1));
            Q_gl = (V_gl(:,1:pd)*Q(:,2:ex+1));
            
            q(:,:,1) = H_gl;
            q(:,:,2) = Q_gl;
            
            if isequal(lim_type, 'minmod')
                [qc,R] = SWE_Limiters.SWEQtoRDG(q,g);
                [R(:,:,1),lim_on_R1] = SWE_Limiters.SlopeLimitCSDG(x_gl_all,R(:,:,1),pd-1,dx,ex,V_gl,inv(V_gl));
                [R(:,:,2),lim_on_R2] = SWE_Limiters.SlopeLimitCSDG(x_gl_all,R(:,:,2),pd-1,dx,ex,V_gl,inv(V_gl));
                q = SWE_Limiters.SWERtoQDG(R,qc,g);
            elseif isequal(lim_type, 'weno')
                [QQ,Xm,Xp] = SWE_Limiters.WENODGWeights(m,inv(V_gl));
                [qc,R] = SWE_Limiters.SWEQtoRDG(q,g);
                [R(:,:,1),lim_on_R1] = SWE_Limiters.WENOlimitDG(x_gl_all,R(:,:,1),pd-1,dx,ex,V_gl,inv(V_gl),QQ,Xm,Xp);
                [R(:,:,2),lim_on_R2] = SWE_Limiters.WENOlimitDG(x_gl_all,R(:,:,2),pd-1,dx,ex,V_gl,inv(V_gl),QQ,Xm,Xp);
                q = SWE_Limiters.SWERtoQDG(R,qc,g);
            end
            
            H_gl = q(:,:,1);
            Q_gl = q(:,:,2);
            
            Q = (inv(V_gl)*Q_gl)';
            H = (inv(V_gl)*H_gl)';
            
            limiter_on(:,1) = lim_on_R1;
            limiter_on(:,2) = lim_on_R2;
        end

    
        function psi = minmod_new(v)
            v=v';
            N = size(v,1);
            m = size(v,2);
            psi = zeros(1,N);
            
            s = sum(sign(v),2)/m;
            ids = find(abs(s)==1);
            
            if (~isempty(ids))
                psi(ids) = s(ids).*min(abs(v(ids,:)),[],2);
            end
        end
            
    
        function [qc,R] = SWEQtoRDG(q,g)
            dim = size(q);
            m = dim(1)-1;
            N = dim(2);
            R = zeros(m+1,N,2);
            qc = zeros(2,N);

            HB = q(:,:,1);
            Q = q(:,:,2);

            qc(1,:)=(1/6)*HB(1,:)+(4/6)*HB(2,:)+(1/6)*HB(3,:);
            qc(2,:)=(1/6)*Q(1,:)+(4/6)*Q(2,:)+(1/6)*Q(3,:);

            for i = 1:N
                [~, iS, ~] = SWE_Limiters.SWEChar([qc(1,i) qc(2,i)],g);
                qh = [HB(:,i) Q(:,i)];
                Ch = (iS*qh')';
                R(:,i,1) = Ch(:,1);
                R(:,i,2) = Ch(:,2);
            end
        end

        function [q] = SWERtoQDG(R,qc,g)
            dim = size(R);
            m = dim(1)-1;
            N = dim(2);
            q = zeros(m+1,N,2);

            for i = 1:N
                [S, ~, ~] = SWE_Limiters.SWEChar([qc(1,i) qc(2,i)],g);
                qh = [R(:,i,1) R(:,i,2)];
                Ch = (S*qh')';
                q(:,i,1) = Ch(:,1);
                q(:,i,2) = Ch(:,2);
            end
        end

        function [S, iS, Lam] = SWEChar(q0, g)
            H0 = q0(1);
            Q0 = q0(2);
            u0 = Q0/H0;
            c0 = sqrt(g*H0);

            S = zeros(2,2);
            Lam = zeros(2,2);

            S(1,1) = 1;
            S(1,2) = 1;
            S(2,1) = u0+c0;
            S(2,2) = u0-c0;

            iS = inv(S);

            Lam(1,1) = u0+c0;
            Lam(2,2) = c0-c0;
        end

        function [ulimit,limiter_on] = SlopeLimitCSDG(x,u,m,h,N,V,iV);
            eps0 = 0.5;
            theta = 1;

            ucell = (1/6)*u(1,:)+(4/6)*u(2,:)+(1/6)*u(3,:);
            ulimit = u;

            ve = [ucell(1) ucell ucell(end)];

            uel = u(1,:);
            uer = u(end,:);
            vj = ucell;
            vjm = ve(1:N);
            vjp = ve(3:N+2);

            vel = vj - SWE_Limiters.minmod_new([vj-uel; vj-vjm; vjp-vj]);
            ver = vj + SWE_Limiters.minmod_new([uer-vj; vj-vjm; vjp-vj]);
            ids = find(abs(vel-uel)>eps0 | abs(ver-uer)>eps0);

            if (~isempty(ids))
                uhl = iV*u(:,ids);
                uhl(3:m+1,:)=0;
                ulin = V*uhl;
                ux = 2/h*(vj(ids)-ulin(1,:));

                x0h = ones(m+1,1)*(x(end,:)+x(1,:))/2; % center of each cell
                ulimit(:,ids) = ones(m+1,1)*vj(ids) + (x(:,ids)-x0h(:,ids)).*(ones(m+1,1)*...
                    SWE_Limiters.minmod_new([ux(1,:); theta*(vjp(ids)-vj(ids))./h; theta*(vj(ids)-vjm(ids))./h]));
            end
            
            limiter_on = zeros(N,1);
            limiter_on(ids) = 1;
        end
        
        function [ulimit, limiter_on] = WENOlimitDG(x,u,m,h,N,V,iV,Q,Xm,Xp);
            eps0 = 1e-8;
            eps1 = 1e-10;
            p = 1;
            gammam1 = 0.001;
            gamma0 = 0.998;
            gammap1 = 0.001;
            
            ucell = (1/6)*u(1,:)+(4/6)*u(2,:)+(1/6)*u(3,:);
            ulimit = u;
            
            ue = [u(:,1) u u(:,end)];
            Pm = Xp'*ue;
            Pp = Xm'*ue;
            Ph = iV*Pm;
            Ph(1,:) = 0;
            Pm = V*Ph;
            Ph = iV*Pp;
            Ph(1,:) = 0;
            Pp = V*Ph;
            
            ve = [ucell(1) ucell ucell(end)];
            
            uel = u(1,:);
            uer = u(end,:);
            vj = ucell;
            vjm = ve(1:N);
            vjp = ve(3:N+2);

            vel = vj - SWE_Limiters.minmod_new([vj-uel; vj-vjm; vjp-vj]);
            ver = vj + SWE_Limiters.minmod_new([uer-vj; vj-vjm; vjp-vj]);
            ids = find(abs(vel-uel)>eps0 | abs(ver-uer)>eps0);

            if (~isempty(ids))
                pm1 = Pm(:,ids) + ones(m+1,1)*vj(ids);
                p0 = u(:,ids);
                pp1 = Pp(:,ids+2) + ones(m+1,1)*vj(ids);
                
                betam1 = diag(pm1'*Q*pm1);
                alpham1 = gammam1./(eps1 + betam1).^(2*p);
                
                beta0 = diag(p0'*Q*p0);
                alpha0 = gamma0./(eps1 + beta0).^(2*p);
                
                betap1 = diag(pp1'*Q*pp1);
                alphap1 = gammap1./(eps1 + betap1).^(2*p);
                
                alphas = alpham1 + alpha0 + alphap1;
                omm1 = alpham1./alphas;
                om0 = alpha0./alphas;
                omp1 = alphap1./alphas;
                
                ulimit(:,ids) = pm1*diag(omm1) + p0*diag(om0) + pp1*diag(omp1);
            end
            
            limiter_on = zeros(N,1);
            limiter_on(ids) = 1;
            
        end
        
        function [Q,Xm,Xp] = WENODGWeights(m,iV)
            Q = zeros(m+1,m+1);
            Pmat = zeros(m+1,m+1);
            Xm = Pmat;
            Xp = Pmat;
            
            [x,w] = SWE_Limiters.LegendreGQ(m);
            Lambda = diag(w);
            
            for i=1:m+1
                Pmat(i,:) = SWE_Limiters.LegendreP(x,i-1)';
                Xm(i,:) = SWE_Limiters.LegendreP(x-2,i-1)';
                Xp(i,:) = SWE_Limiters.LegendreP(x+2,i-1)';
            end
            
            for l=1:m
                H = zeros(m+2-l,m+2-l);
                H(1,1) = 1/sqrt((2*l+1)*(2*l-1));
                H(m+2-l,m+2-l) = 1/(sqrt(2*(m+2)+1)*sqrt(2*(m+2)-1));
                for i=2:m-l+1
                    Ah = 1/(sqrt(2*(l-1+i)+1)*sqrt(2*(l-1+i)-1));
                    H(i,i) = Ah;
                    H(i+1,i-1) = -Ah;
                end
                
                Ph1 = H\Pmat(l:m+1,:);
                Pmat(1:l,:)=0;
                Pmat(l+1:m+1,:) = Ph1(1:m-l+1,:);
                
                Qh = Pmat*Lambda*Pmat';
                Q = Q + 2^(2*l-1)*Qh;
            end
            
            Q = iV'*Q*iV;
            
            Xp = iV'*Xp;
            Xm = iV'*Xm;
        end
        
        function P = LegendreP(x,m)
            xp = x;
            dims = size(xp);
            if dims(2) == 1
                xp = xp';
            end
            
            PL = zeros(m+1,length(xp));
            PL(1,:) = sqrt(1.0/2.0);
            if m == 0
                P = PL';
                return;
            end
            PL(2,:) = sqrt(3.0/2.0)*xp;
            if m == 1
                P = PL(m+1,:)';
                return;
            end
            
            aold = sqrt(1.0/3.0);
            for i=1:m-1
                anew = 2/(2*i+2)*sqrt((i+1)*(i+1)*(i+1)*(i+1)/(2*i+1)/(2*i+3));
                PL(i+2,:) = 1/anew*(-aold*PL(i,:) + xp.*PL(i+1,:));
                aold=anew;
            end
            P = PL(m+1,:)';
        end
        
        function [x,w] = LegendreGQ(m)
            if m == 0
                x(1) = 0;
                w(1) = 2;
                return;
            end
            
            J = zeros(m+1);
            h1 = 2*(0:m);
            J = diag(2./(h1(1:m)+2).*...
                sqrt((1:m).*((1:m)).*((1:m)).*((1:m))./(h1(1:m)+1)./(h1(1:m)+3)),1);
            J(1,1) = 0;
            J = J+J';
            
            [V,D] = eig(J);
            x = diag(D);
            w = 2*(V(1,:)').^2;
        end
    end
    
end