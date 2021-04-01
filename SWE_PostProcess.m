classdef SWE_PostProcess
    methods(Static)
        
        % Plot function
        function plot_pIteratations(m_iter,plot_times,P,pnew,data)
            
            set(0, 'defaultAxesTickLabelInterpreter','latex'); 
            set(0, 'defaultLegendInterpreter','latex');
            set(0, 'defaultTextInterpreter', 'latex');
            
            choice = data.choice;
            filename_str = data.filename_str;
            
            Pfunc = data.p_true;
            Ptrue = Pfunc(data.t_steps_meas);
            
            % Plotting colors
            grey = (1/255)*[166, 166, 166];
            light_green = (1/255)*[152, 230, 152];
            light_blue = (1/255)*[102, 204, 255];
            light_purple = (1/255)*[204, 153, 255];
            brown = (1/255)*[153, 102, 51];
            red = (1/255)*[255, 51, 51];
            orange = (1/255)*[255, 133, 51];
            yellow = (1/255)*[255, 204, 0];
            green = (1/255)*[51, 153, 102];
            blue  = (1/255)*[0, 102, 204];
            purple = (1/255)*[128, 0, 128];
            allcolors = [light_green; light_blue; light_purple; brown; red; orange; yellow; green; blue; purple;...
                light_green; light_blue; light_purple; brown; red; orange; yellow; green; blue; purple];
            
            % Font & Marker Sizes
            title_fs = 45;
            subtitle_fs = 40;
            xlabel_fs = 35;
            legend_fs = 30;
            sublegend_fs = 20;
            ticks_fs = 30;
            marker_s1 = 500;
            marker_s2 = 250;
            
            % Plot
            fig1 = figure(1);
            fig2 = figure(2);
            
           if mod(m_iter,plot_times) == 0
                set(0, 'CurrentFigure', fig1)
                scatter(([m_iter-plot_times+1:m_iter]), (data.iter_err1(m_iter-plot_times+1:m_iter)), marker_s2, green, 'LineWidth', 3)
                hold on;
                if m_iter == data.num_of_iters
                    scatter((data.best_iter1), (data.best_err1), marker_s1, red, '+', 'LineWidth', 3)
                    hold on;
                end

                xl1 = xlabel('Iteration');
                xl1.FontSize = xlabel_fs;
                yl1 = ylabel('Iteration Error');
                yl1.FontSize = xlabel_fs;
                ax = gca;
                ax.XAxis.FontSize = ticks_fs;
                ax.YAxis.FontSize = ticks_fs; 
                
                h = zeros(2, 1);
                h(1) = scatter(NaN,NaN, marker_s2, green, 'LineWidth', 3);
                h(2) = scatter(NaN,NaN, marker_s1, red, '+', 'LineWidth', 3);
                lg = legend(h, 'Iteration Error', 'Best Iteration Error');
                lg.Location = 'northeast';
                lg.FontSize = xlabel_fs;
                legend boxoff
                set(gca,'yscale','log','xscale','log')
                drawnow()
           end
            

            
            set(0, 'CurrentFigure', fig2)
            subplot(2,3,[1 2 4 5])
            if m_iter == 1
                plot(data.t_steps_meas, Ptrue, 'k', 'LineWidth', 3, 'DisplayName','True')
                hold on;
                plot(data.t_steps_forward, P, '--', 'Color', grey,'LineWidth', 3, 'DisplayName','Noisy Initial Guess')
                hold on;
            end
            if mod(m_iter,plot_times) == 0
                if m_iter <= 10*plot_times
                    plot(data.t_steps_forward,pnew, ':', 'Color', allcolors(m_iter/plot_times,:), 'LineWidth', 3, 'DisplayName',['Iteration ', num2str(m_iter)])
                else
                    plot(data.t_steps_forward,pnew, '-.', 'Color', allcolors(m_iter/plot_times,:), 'LineWidth', 3, 'DisplayName',['Iteration ', num2str(m_iter)])
                end
                hold on;
                t2 = title('P');
                t2.FontSize = title_fs;
                xlabel('time', 'FontSize', xlabel_fs)
%                 legend('Location','eastoutside', 'FontSize', legend_fs);
                ax = gca;
                ax.XAxis.FontSize = ticks_fs;
                ax.YAxis.FontSize = ticks_fs;
                drawnow()
            end
            
            subplot(2,3,[3 6])
            if m_iter == 1
                plot(NaN, NaN,'k', 'LineWidth', 2, 'DisplayName',['True'])
                hold on;
            end
            if mod(m_iter,plot_times) == 0
                if m_iter <= 10*plot_times
                    plot(NaN, NaN, ':', 'Color', allcolors(m_iter/plot_times,:), 'LineWidth', 2, 'DisplayName',['Iteration ', num2str(m_iter)])
                else
                    plot(NaN, NaN, '-.', 'Color', allcolors(m_iter/plot_times,:), 'LineWidth', 2, 'DisplayName',['Iteration ', num2str(m_iter)])
                end
                hold on;
                set(gca,'xtick',[], 'ytick', [])
                drawnow()
            end
                
            
            if mod(m_iter,plot_times) == 0 || m_iter == 1
                
                set(0, 'CurrentFigure', fig1)
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                                
                set(0, 'CurrentFigure', fig2)
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
            end
            
            % Save
            if mod(m_iter,plot_times) == 0
                plot_filename1jpg = sprintf('%s/SWEInversePlot_IterationErrors_%s.jpg',data.front_path,filename_str);
                plot_filename1fig = sprintf('%s/SWEInversePlot_IterationErrors_%s.fig',data.front_path,filename_str);
                saveas(fig1,plot_filename1jpg)
                saveas(fig1,plot_filename1fig)
                
                plot_filename2jpg = sprintf('%s/SWEInversePlot_PIterations_%s.jpg',data.front_path,filename_str);
                plot_filename2fig = sprintf('%s/SWEInversePlot_PIterations_%s.fig',data.front_path,filename_str);
                saveas(fig2,plot_filename2jpg)
                saveas(fig2,plot_filename2fig)
            end
        end
        

        % Error function
        function errors(ex_all, HH, QQ, delx, T, pd, choice)


            %% ------------------------------------------------------------------------
            %     Initialize Vectors for Accuracy Test
            % -------------------------------------------------------------------------

            iterations = length(ex_all);
            initialex = ex_all(1);

            % For Storing L^1 & L^2 Calculations
            L1ErrorH  = zeros(iterations-1,1);
            [L1QuotientH, L1AccuracyH, L1ErrorQ, L1QuotientQ, L1AccuracyQ,...
                L2ErrorH, L2QuotientH, ...
                L2AccuracyH, L2ErrorQ, L2QuotientQ, L2AccuracyQ, LinfErrorH, LinfQuotientH, LinfAccuracyH, ...
                LinfErrorQ, LinfQuotientQ, ex_size] = deal(L1ErrorH);


            %% ------------------------------------------------------------------------
            %     Calculate Error
            % -------------------------------------------------------------------------

            pwrs = 0:pd-1;
            % Quadrature points in reference cell
            num_nodes = 5;
            [xi, wi] = lgwt(num_nodes,-0.5,0.5);
            xi = sort(xi);


            for r = 2:iterations
                ex = initialex*2^(r-2);
                ex_size(r-1) = ex;
                ratio = 2;
                H_ex_quad = zeros(2*ex,1);
                Q_ex_quad = zeros(2*ex,1);

                % Quadrature points and weights
                [xi_2ex,wi_2ex] = lgwt(num_nodes, -0.5, 0.5);
                [xi_ex_neg,~] = lgwt(num_nodes, -0.5, 0);
                [xi_ex_pos,~] =  lgwt(num_nodes, 0, 0.5);

                H_ex = squeeze(HH(r-1,:,1:ex));
                H_2ex = squeeze(HH(r,:,1:ratio*ex));                    

                Q_ex = squeeze(QQ(r-1,:,1:ex));
                Q_2ex = squeeze(QQ(r,:,1:ratio*ex));

                if pd == 1
                    H_ex_quad_neg = ((xi_ex_neg.^pwrs)*H_ex');
                    Q_ex_quad_neg = ((xi_ex_neg.^pwrs)*Q_ex');

                    H_ex_quad_pos = ((xi_ex_pos.^pwrs)*H_ex');
                    Q_ex_quad_pos = ((xi_ex_pos.^pwrs)*Q_ex');
                else
                    H_ex_quad_neg = ((xi_ex_neg.^pwrs)*H_ex);
                    Q_ex_quad_neg = ((xi_ex_neg.^pwrs)*Q_ex);

                    H_ex_quad_pos = ((xi_ex_pos.^pwrs)*H_ex);
                    Q_ex_quad_pos = ((xi_ex_pos.^pwrs)*Q_ex);
                end

                H_ex_quad = H_ex_quad_pos(:,[1;1]*(1:size(H_ex_quad_pos,2)));
                Q_ex_quad = Q_ex_quad_pos(:,[1;1]*(1:size(Q_ex_quad_pos,2)));

                H_ex_quad(:,1:2:end) = H_ex_quad_neg;
                Q_ex_quad(:,1:2:end) = Q_ex_quad_neg;

                if pd == 1
                    H_2ex_quad = ((xi_2ex.^pwrs)*H_2ex');
                    Q_2ex_quad = ((xi_2ex.^pwrs)*Q_2ex');
                else
                    H_2ex_quad = ((xi_2ex.^pwrs)*H_2ex);
                    Q_2ex_quad = ((xi_2ex.^pwrs)*Q_2ex);
                end

                HtempL1(r-1,1:2*ex) = sum(abs(H_2ex_quad-H_ex_quad).*wi_2ex);
                QtempL1(r-1,1:2*ex) = sum(abs(Q_2ex_quad-Q_ex_quad).*wi_2ex);

                HtempL2(r-1,1:2*ex) = sum(((H_2ex_quad-H_ex_quad).^2).*wi_2ex);
                QtempL2(r-1,1:2*ex) = sum(((Q_2ex_quad-Q_ex_quad).^2).*wi_2ex);

                L1ErrorH(r-1) = delx(r)*sum(HtempL1(r-1,:));
                L1ErrorQ(r-1) = delx(r)*sum(QtempL1(r-1,:));

                L2ErrorH(r-1) = sqrt(delx(r)*sum(HtempL2(r-1,:)));
                L2ErrorQ(r-1) = sqrt(delx(r)*sum(QtempL2(r-1,:)));

            end
                


            for r = 1:iterations-2 
                rr = r+1;

                % L^1 Accuracy
                L1QuotientH(rr) = L1ErrorH(r)/L1ErrorH(r+1);
                L1AccuracyH(rr) = log2(L1QuotientH(rr));
                L1QuotientQ(rr) = L1ErrorQ(r)/L1ErrorQ(r+1);
                L1AccuracyQ(rr) = log2(L1QuotientQ(rr));
                
                % L^2 Accuracy
                L2QuotientH(rr) = L2ErrorH(r)/L2ErrorH(r+1);
                L2AccuracyH(rr) = log2(L2QuotientH(rr));
                L2QuotientQ(rr) = L2ErrorQ(r)/L2ErrorQ(r+1);
                L2AccuracyQ(rr) = log2(L2QuotientQ(rr));
                
            end

            %% ------------------------------------------------------------------------
            %     Display Error Matrices
            % -------------------------------------------------------------------------


            fprintf(['L^1 Error Table at T = ', num2str(T)])
            ErrorTableL1 = table(ex_size,L1ErrorH,L1AccuracyH,...
                L1ErrorQ,L1AccuracyQ);
            ErrorTableL1.Properties.VariableNames = ...
                 {'N', 'L1_Error_H','L1_Order_H','L1_Error_Q','L1_Order_Q'};
            fprintf('\n')
            disp(ErrorTableL1)

            fprintf(['L^2 Error Table at T = ', num2str(T)])
            ErrorTableL2 = table(ex_size,L2ErrorH,L2AccuracyH,...
                L2ErrorQ,L2AccuracyQ);
            ErrorTableL2.Properties.VariableNames = ...
                 {'N', 'L2_Error_H','L2_Order_H','L2_Error_Q','L2_Order_Q'};
            fprintf('\n')
            disp(ErrorTableL2)


        end

    end
end
