function [best_error, best_reg] = SWE_Results(k,nt,func_val,meshes,hyparams,accuracy_test)

    close all
    
    % Basic Info
    choice = str2double(sprintf('%d',func_val));
    data = SWE_Data.data_params(k,nt,func_val,meshes,hyparams,choice,accuracy_test);
    filename_str = data.filename_str;  
    
    % Use Latex formatting
    set(0, 'defaultAxesTickLabelInterpreter','latex'); 
    set(0, 'defaultLegendInterpreter','latex');
    set(0, 'defaultTextInterpreter', 'latex');

    % Plotting Colors
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
    
    % Plotting Line Styles
    lstyles = ["--", ":", "-.", "-"];

    % Font & Marker Sizes
    title_fs = 85;
    subtitle_fs = 50;
    xlabel_fs = 75;
    legend_fs = 75;
    sublegend_fs = 35;
    ticks_fs = 70;
    marker_s1 = 800;
    marker_s2 = 250;
    
    % Read in Files
    file_errors = sprintf('%s/SWE_IterationErrors_%s.csv',data.front_path,data.filename_str);
    file_besterror = sprintf('%s/SWE_BestError_%s.csv',data.front_path,data.filename_str);
    file_bestreg = sprintf('%s/SWE_BestReg_%s.csv',data.front_path,data.filename_str);
    file_errorsreg = sprintf('%s/SWE_IterationErrorsReg_%s.csv',data.front_path,data.filename_str);
    file_allp = sprintf('%s/SWE_PData_%s.csv',data.front_path,data.filename_str);
    file_bestp = sprintf('%s/SWE_BestP_%s.csv',data.front_path,data.filename_str);
    file_bestH = sprintf('%s/SWE_BestH_%s.csv',data.front_path,data.filename_str);
    file_bestQ = sprintf('%s/SWE_BestQ_%s.csv',data.front_path,data.filename_str);
    file_bestB = sprintf('%s/SWE_BestB_%s.csv',data.front_path,data.filename_str);
    file_trueH = sprintf('%s/SWE_Meas_%s_H.csv',data.front_path,data.filename_str);
    file_trueQ = sprintf('%s/SWE_Meas_%s_Q.csv',data.front_path,data.filename_str);
    file_trueB = sprintf('%s/SWE_Meas_%s_B.csv',data.front_path,data.filename_str);
    file_tsteps_meas = sprintf('%s/SWE_tsteps_meas_%s.csv',data.front_path,data.filename_str);
    file_tsteps_forward = sprintf('%s/SWE_ForwardTimeSteps_%s.csv',data.front_path,data.filename_str);
    file_pguess = sprintf('%s/SWE_PGuessNoisy_%s.csv',data.front_path,data.filename_str);
    
    errors = csvread(file_errors);
    best_err = csvread(file_besterror);
    best_reg = csvread(file_bestreg);
    if isfile(file_errorsreg)
        errors_reg = csvread(file_errorsreg);
    end
    all_p = csvread(file_allp);
    best_p = csvread(file_bestp);
    best_h = csvread(file_bestH);
    best_q = csvread(file_bestQ);
    best_b = csvread(file_bestB);
    true_h = csvread(file_trueH);
    true_q = csvread(file_trueQ);
    true_b = csvread(file_trueB);
    t_steps_meas = csvread(file_tsteps_meas);
    t_steps_forward = csvread(file_tsteps_forward);
    p_guess = csvread(file_pguess);
    p_true = data.p_true(t_steps_meas);
    best_error = errors(best_err);
    
    % Adjust/Reshape Data
    nt_meas = length(t_steps_meas);
    nt_forward = length(t_steps_forward);
    ratio = ceil(nt_meas/nt_forward);
    shift = nt_forward - mod(nt_meas,nt_forward);
    save_pts_forward = [floor((1/4)*nt_forward), floor((2/4)*nt_forward), floor((3/4)*nt_forward), nt_forward];
    save_pts_meas = save_pts_forward*ratio - shift;
    best_h = reshape(best_h,size(best_h,1),k(2)+1,nt_forward);
    best_q = reshape(best_q,size(best_q,1),k(2)+1,nt_forward);
    best_b = reshape(best_b,size(best_b,1),k(2)+1,nt_forward);
    best_h = best_h(:,:,save_pts_forward);
    best_q = best_q(:,:,save_pts_forward);
    best_b = best_b(:,:,save_pts_forward);
    
    true_h = reshape(true_h,size(true_h,1),k(1)+1,nt_meas);
    true_q = reshape(true_q,size(true_q,1),k(1)+1,nt_meas);
    true_b = reshape(true_b,size(true_b,1),k(1)+1,nt_meas);
    true_h = true_h(:,:,save_pts_meas);
    true_q = true_q(:,:,save_pts_meas);
    true_b = true_b(:,:,save_pts_meas);
        
    % Get Cell Averages
    best_h_ca = zeros(length(save_pts_forward),size(best_h,1));
    [best_q_ca, best_b_ca] = deal(best_h_ca);
    if k(2) == 0 || k(2) == 1
        best_h_ca = squeeze(best_h(:,1,:));
        best_q_ca = squeeze(best_q(:,1,:));
        best_b_ca = squeeze(best_b(:,1,:));
    elseif k(2) == 2 || k(2) == 3
        for i = 1:length(save_pts_forward)
            best_h_ca(i,:) = squeeze(best_h(:,1,i) + 1/12*best_h(:,3,i))';
            best_q_ca(i,:) = squeeze(best_q(:,1,i) + 1/12*best_q(:,3,i))';
            best_b_ca(i,:) = squeeze(best_b(:,1,i) + 1/12*best_b(:,3,i))';
        end
    end
    true_h_ca = zeros(length(save_pts_meas),size(true_h,1));
    [true_q_ca, true_b_ca] = deal(true_h_ca);
    if k(1) == 0 || k(1) == 1
        true_h_ca = squeeze(true_h(:,1,:));
        true_q_ca = squeeze(true_q(:,1,:));
        true_b_ca = squeeze(true_b(:,1,:));
    elseif k(1) == 2 || k(1) == 3
        for i = 1:length(save_pts_meas)
            true_h_ca(i,:) = squeeze(true_h(:,1,i) + 1/12*true_h(:,3,i))';
            true_q_ca(i,:) = squeeze(true_q(:,1,i) + 1/12*true_q(:,3,i))';
            true_b_ca(i,:) = squeeze(true_b(:,1,i) + 1/12*true_b(:,3,i))';
        end
    end
    
    % Plot initializations
    fig1 = figure(1);
    fig2 = figure(2);
    fig3 = figure(3);
    fig4 = figure(4);
    fig5 = figure(5);
    fig6 = figure(6);

    % PLOT 1: Error Plots (w/o Regularization) - LogLog Scale
    set(0, 'CurrentFigure', fig1)
    scatter(1:length(errors), errors, marker_s2, green, 'LineWidth', 3)
    hold on;
    scatter(best_err, errors(best_err), marker_s1, red, '+', 'LineWidth', 3)
    
    xlabel('Iteration', 'FontSize', xlabel_fs);
    ylabel('Iteration Error', 'FontSize', xlabel_fs);
    ax = gca;
    ax.XAxis.FontSize = ticks_fs;
    ax.YAxis.FontSize = ticks_fs; 
    [~, objh] = legend({'Iteration Error', ['Smallest Error: Iteration ', num2str(best_err)]}, 'Location', 'northoutside', 'Fontsize', legend_fs);
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 20, 'LineWidth', 3); %// set marker size as desired
    objhl = findobj(objh, 'type', 'patch'); % objects of legend of type patch
    set(objhl, 'Markersize', 20, 'LineWidth', 3); % set marker size as desired
    legend boxoff
    set(gca,'yscale','log','xscale','log')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    plot_filename1jpg = sprintf('%s/SWEPlot_LogLogErrors_%s.jpg',data.front_path,filename_str);
    plot_filename1fig = sprintf('%s/SWEPlot_LogLogErrors_%s.fig',data.front_path,filename_str);
    saveas(fig1,plot_filename1jpg)
    saveas(fig1,plot_filename1fig)
    
    if isfile(file_errorsreg)
        % PLOT 2: Error Plots (w/ Regularization) 
        set(0, 'CurrentFigure', fig2)

        % START USER CHOICES --------------------------------------------------
        fig_type = 4;

        if fig_type == 1
            % Case 1: Global picture with log(n) on x-axis and no fit line
            start_val = 1;
            end_val = 1000;
            x_type = 'log';
            y_type = 'err_last';
            y_scale = 'log';
            overlay_line = 0;
            overlay_start = start_val;
        elseif fig_type == 2
            % Case 2: Global picture with n on x-axis and no fit line
            start_val = 1;
            end_val = 1000;
            x_type = 'standard';
            y_type = 'err_last';
            y_scale = 'log';
            overlay_line = 0;
            overlay_start = start_val;
        elseif fig_type == 3
            % Case 3: Front behavior with log(n) on x-axis and a fit line
            start_val = 1;
    %         end_val = 50;
    %         end_val = 40;
            end_val = 20;
    %         end_val = 15;
            x_type = 'log';
            y_type = 'err_last';
            y_scale = 'log';
            overlay_line = 1;
            overlay_start = start_val;
        elseif fig_type == 4
            % Case 3: Tail behavior with n on x-axis and a fit line
            start_val = 50;
            end_val = 1000;
            x_type = 'standard';
            y_type = 'err_last';
            y_scale = 'log';
            overlay_line = 1;
            overlay_start = start_val;
        elseif fig_type == 5
            % Customize with Other options to select:
            % CHOICES 1-2: Range of iterations to plot
            start_val = 1;
        %     start_val = 50;
        %     start_val = 600;
        %     end_val = 1000;
        %     end_val = 1500;
            end_val = 50;

            % CHOICE 3: Choose standard or log for x-axis
        %     x_type = 'standard';
            x_type = 'log';

            % CHOICE 4: Choose valuse for plotting on the y_axis
        %     y_type = 'err_only';
            y_type = 'err_last';
        %     y_type = 'err_min';

            % CHOICE 5: Choose standard or log for y-axis
        %     y_scale = 'standard';
            y_scale = 'log';

            % CHOICES 6-7 Choose if want a best fit line and where to start it
        %     overlay_line = 0;
            overlay_line = 1;

        %     overlay_start = 1;
            overlay_start = start_val;
        end
        % END USER CHOICES ----------------------------------------------------

        % Determine x values based on input
        if strcmp(x_type, 'standard') 
            x_plot = (start_val:end_val);
            x_label = 'Iteration n';
        elseif strcmp(x_type, 'log') 
            x_plot = log(start_val:end_val);
            x_label = 'log(n)';
        end

        % Determine y values based on inputs
        if strcmp(y_scale, 'standard')
            if strcmp(y_type, 'err_only') 
                y_plot = errors_reg(start_val:end_val);
                y_label = 'E(n)';
            elseif strcmp(y_type, 'err_last') 
                y_plot = errors_reg(start_val:end_val) - errors_reg(end);
                y_label = '$E(n) - E(*)$';
            elseif strcmp(y_type, 'err_min') 
                y_plot = errors_reg(start_val:end_val) - min(errors_reg);
                y_label = 'E(n) - min{E(n)}';
            end
        elseif strcmp(y_scale, 'log')
            if strcmp(y_type, 'err_only') 
                y_plot = log(errors_reg(start_val:end_val));
                y_label = 'log(E(n))';
            elseif strcmp(y_type, 'err_last') 
                y_plot = log(errors_reg(start_val:end_val) - errors_reg(end));
                y_label = 'log(\textbf{J}(n) - \textbf{J}(*))';
            elseif strcmp(y_type, 'err_min') 
                y_plot = log(errors_reg(start_val:end_val) - min(errors_reg));
                y_label = 'log(E(n) - min{E(n)})';
            end
        end

        % Values for linear fit line if selected
        if overlay_line == 1
            c = polyfit(x_plot, y_plot, 1);
            if strcmp(x_type, 'standard') 
                x_est = (overlay_start:end_val);
            elseif strcmp(x_type, 'log') 
                x_est = log(overlay_start:end_val);
            end
            y_est = polyval(c, x_est);
            disp(c)
        end

        % Create plot
        scatter(x_plot, y_plot, marker_s2, green, 'LineWidth', 3)
        hold on;
        if overlay_line == 1
            plot(x_est, y_est, 'Color', blue, 'LineWidth', 3)
        end

        % Format the plots
        xlabel(x_label, 'FontSize', xlabel_fs);
        ylabel(y_label, 'FontSize', xlabel_fs);
        ax = gca;
        ax.XAxis.FontSize = ticks_fs;
        ax.YAxis.FontSize = ticks_fs; 
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        plot_filename2jpg = sprintf('%s/SWEPlot_ErrorsRegCase%d_%s.jpg',data.front_path,fig_type,filename_str);
        plot_filename2fig = sprintf('%s/SWEPlot_ErrorsRegCase%d_%s.fig',data.front_path,fig_type,filename_str);
        saveas(fig2,plot_filename2jpg)
        saveas(fig2,plot_filename2fig)
    end
    
    % PLOT 3: Best Choice of P
    set(0, 'CurrentFigure', fig3)
    plot(t_steps_meas, p_true, 'k', 'LineWidth', 8)
    hold on;
    plot(t_steps_forward, best_p, '--', 'Color', green,'LineWidth', 8)
    xlabel('time', 'FontSize', xlabel_fs)
    xlim([0,data.T])
    ylim([min(min(best_p),min(p_true)),1.2*max(max(p_true),max(best_p))])
    ax = gca;
    ax.XAxis.FontSize = ticks_fs;
    ax.YAxis.FontSize = ticks_fs;
    lg = legend('$p_{true}(t)$', ['Best Iteration: $p^{', num2str(best_err), '}(t)$'], 'Location', 'northeast');
    lg.Location = 'northoutside';
    lg.FontSize = legend_fs;
    legend boxoff;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    plot_filename3jpg = sprintf('%s/SWEPlot_BestP_%s.jpg',data.front_path,filename_str);
    plot_filename3fig = sprintf('%s/SWEPlot_BestP_%s.fig',data.front_path,filename_str);
    saveas(fig3,plot_filename3jpg)
    saveas(fig3,plot_filename3fig)
    
    % PLOT 4: True H vs Best H
    set(0, 'CurrentFigure', fig4)
    lg_text_meas_h = ["$(h+b)_{true}(x,\frac{T}{3})$", "$(h+b)_{true}(x,\frac{2T}{3})$", "$(h+b)_{true}(x,T)$"];
    lg_text_meas_b = ["$b_{true}(x,\frac{T}{3})$", "$b_{true}(x,\frac{2T}{3})$", "$b_{true}(x,T)$"];
    lg_text_best_h = ["$(h+b)_{best}(x,\frac{T}{3})$", "$(h+b)_{best}(x,\frac{2T}{3})$", "$(h+b)_{best}(x,T)$"];
    lg_text_best_b = ["$b_{best}(x,\frac{T}{3})$", "$b_{best}(x,\frac{2T}{3})$", "$b_{best}(x,T)$"];
    title_text_time = ["$t = \frac{T}{4}$", "$t = \frac{T}{2}$", "$t = \frac{3T}{4}$", "$t=T$"];
    
    fig4_t = tiledlayout(1,length(save_pts_meas)+1,'Padding','Compact');
    
    for i = 1:length(save_pts_meas)
        nexttile
        plot(data.meas_mesh_center(1:end-1),true_h_ca(i,:)+true_b_ca(i,:),lstyles(4),'Color', 'k', 'LineWidth', 3, 'DisplayName', '$(h+b)_{meas}(x,t)$')
        hold on;
        plot(data.meas_mesh_center(1:end-1),true_b_ca(i,:),lstyles(4),'Color', blue, 'LineWidth', 3, 'DisplayName', '$(b)_{meas}(x,t)$')
        hold on;
    
        plot(data.forward_mesh_center(1:end-1),best_h_ca(i,:)+best_b_ca(i,:),lstyles(2), 'Color', green, 'LineWidth', 3, 'DisplayName', '$(h+b)_{best}(x,t)$')
        hold on;
        plot(data.forward_mesh_center(1:end-1),best_b_ca(i,:),lstyles(2), 'Color', red, 'LineWidth', 3, 'DisplayName', '$(b)_{best}(x,t)$')
        hold on;
        
        xlabel('x', 'FontSize', sublegend_fs)
        if i == 1
            ylabel({' ';'Water Surface Height,'; 'Bottom Topography'}, 'FontSize', sublegend_fs)
        else
            ylabel({' ';' '}, 'FontSize', sublegend_fs)
        end
        title(title_text_time(i), 'FontSize', sublegend_fs)
        ylim([0,15])
        ax = gca;
        ax.XAxis.FontSize = sublegend_fs;
        ax.YAxis.FontSize = sublegend_fs;
        axis square;
    end
    nexttile
    plot(NaN,NaN,lstyles(4),'Color', 'k', 'LineWidth', 3, 'DisplayName', '$(h+b)_{meas}(x,t)$')
    hold on;
    plot(NaN,NaN,lstyles(4),'Color', blue, 'LineWidth', 3, 'DisplayName', '$(b)_{meas}(x,t)$')
    hold on;

    plot(NaN,NaN,lstyles(2), 'Color', green, 'LineWidth', 3, 'DisplayName', '$(h+b)_{best}(x,t)$')
    hold on;
    plot(NaN,NaN,lstyles(2), 'Color', red, 'LineWidth', 3, 'DisplayName', '$(b)_{best}(x,t)$')
    hold on;
    Lgnd = legend('show');
    Lgnd.Location = 'west';
    Lgnd.FontSize = 30;
    legend boxoff;
    axis off;

    set(gcf,'units','normalized','outerposition',[0 0 1 0.5])
    plot_filename4jpg = sprintf('%s/SWEPlot_ForwardCompareH_%s.jpg',data.front_path,filename_str);
    plot_filename4fig = sprintf('%s/SWEPlot_ForwardCompareH_%s.fig',data.front_path,filename_str);
    saveas(fig4,plot_filename4jpg)
    saveas(fig4,plot_filename4fig)
    
    % PLOT 5: True Q vs Best Q
    set(0, 'CurrentFigure', fig5)
    lg_text_meas_hu = ["$(hu)_{true}(x,\frac{T}{3})$", "$(hu)_{true}(x,\frac{2T}{3})$", "$(hu)_{true}(x,T)$"];
    lg_text_best_hu = ["$(hu)_{best}(x,\frac{T}{3})$", "$(hu)_{best}(x,\frac{2T}{3})$", "$(hu)_{best}(x,T)$"];
    title_text_time = ["$t = \frac{T}{4}$", "$t = \frac{T}{2}$", "$t = \frac{3T}{4}$", "$t=T$"];
    
    fig5_t = tiledlayout(1,length(save_pts_meas)+1,'Padding','Compact');
    
    for i = 1:length(save_pts_meas)
        nexttile
        plot(data.meas_mesh_center(1:end-1),true_q_ca(i,:),lstyles(4),'Color', 'k', 'LineWidth', 3, 'DisplayName', '$(hu)_{meas}(x,t)   $   ')
        hold on;
    
        plot(data.forward_mesh_center(1:end-1),best_q_ca(i,:),lstyles(2), 'Color', green, 'LineWidth', 3, 'DisplayName', '$(hu)_{best}(x,t)   $   ')
        hold on;
        
        xlabel('x', 'FontSize', sublegend_fs)
        if i == 1
            ylabel({' ';'Water Discharge'}, 'FontSize', sublegend_fs)
        else
            ylabel({' ';' '}, 'FontSize', sublegend_fs)
        end
        title(title_text_time(i), 'FontSize', sublegend_fs)
        ylim([-15,20])
        ax = gca;
        ax.XAxis.FontSize = sublegend_fs;
        ax.YAxis.FontSize = sublegend_fs;
        axis square;
    end
    nexttile
    plot(NaN, NaN, lstyles(4),'Color', 'k', 'LineWidth', 3, 'DisplayName', '$(hu)_{meas}(x,t)$')
    hold on;
    plot(NaN, NaN, lstyles(2), 'Color', green, 'LineWidth', 3, 'DisplayName', '$(hu)_{best}(x,t)$')
    Lgnd = legend('show');
    Lgnd.Location = 'west';
    Lgnd.FontSize = 30;
    legend boxoff;
    axis off;
    
    set(gcf,'units','normalized','outerposition',[0 0 1 0.5])
    plot_filename5jpg = sprintf('%s/SWEPlot_ForwardCompareQ_%s.jpg',data.front_path,filename_str);
    plot_filename5fig = sprintf('%s/SWEPlot_ForwardCompareQ_%s.fig',data.front_path,filename_str);
    saveas(fig5,plot_filename5jpg)
    saveas(fig5,plot_filename5fig)
    
    if data.num_of_iters > 500
        % Iterations of P
        if data.num_of_iters == 500
            num_plot_vec = [10,50,100,250,500];
        elseif data.num_of_iters == 750
            num_plot_vec = [10,50,100,250,500,750];
        elseif data.num_of_iters == 1000
            num_plot_vec = [10,50,100,250,500,750,1000];
        elseif data.num_of_iters == 1500
            num_plot_vec = [10,50,100,250,500,750,1000,1500];
        end

        set(0, 'CurrentFigure', fig6)
        plot(t_steps_meas, p_true, 'k', 'LineWidth', 8, 'DisplayName','$p_{true}(t)$')
        hold on;
        plot(t_steps_forward, p_guess, '--', 'Color', grey,'LineWidth', 6, 'DisplayName','$p^0(t)$')
        hold on;

        for i = 1:length(num_plot_vec)
            plot(t_steps_forward, all_p(num_plot_vec(i),:), ':', 'Color', allcolors(i,:), 'LineWidth', 8, 'DisplayName',sprintf('$p^{%d}(t)$',num_plot_vec(i)))
            hold on;
        end
        xlabel('time', 'FontSize', xlabel_fs)
        xlim([0,data.T])
        ax = gca;
        ax.XAxis.FontSize = ticks_fs;
        ax.YAxis.FontSize = ticks_fs;
        lg = legend('Location', 'eastoutside');
        lg.FontSize = legend_fs;
        legend boxoff;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        plot_filename6jpg = sprintf('%s/SWEPlot_Iters_%s.jpg',data.front_path,filename_str);
        plot_filename6fig = sprintf('%s/SWEPlot_Iters_%s.fig',data.front_path,filename_str);
        saveas(fig6,plot_filename6jpg)
        saveas(fig6,plot_filename6fig)
    end
    
end
    
