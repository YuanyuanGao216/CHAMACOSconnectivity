classdef Connectivity < handle
    properties
        DAP_table
        task_flag = {}
        n_task = 0;

        T_ROI_position

        independent_var
        formula

        model_table = table;
    end
    
    methods
        function obj = Connectivity()
        end
        function set_DAP_table(obj,varargin)
            if nargin == 1
                DAP_file = fullfile('..','data','chamacos_fNIRS_variables.xlsx');
            elseif nargin == 2
                DAP_file = varargin{1};
            end
            DAP_table = readtable(DAP_file);
            DAP_table.sex = categorical(DAP_table.sex);
            DAP_table.educcat_mom = categorical(DAP_table.educcat_mom);
            DAP_table.usyr3 = categorical(DAP_table.usyr3);
            obj.DAP_table = DAP_table;
        end
        function set_T_ROI_position(obj)
            obj.T_ROI_position = readtable(fullfile('..','data','T_ROI_position.csv'));
        end

        function stratify(obj, index)
            obj.DAP_table = obj.DAP_table(index,:);
        end

        function add_task(obj,connectivity_matrix,task_flag)
            obj.task_flag{end+1} = task_flag;
            obj.n_task = obj.n_task + 1;

            for i = 1:size(connectivity_matrix,1)
                for j = 1:size(connectivity_matrix,2)
                    if i >= j
                        continue
                    end
                    fc = connectivity_matrix(i,j,:);
                    fc = reshape(fc,[],1);
                    name = [task_flag, '_', num2str(i), '_', num2str(j)];
                    obj.DAP_table.(name) = fc;
                end
            end
        end

        function get_table_names(obj)
            for i = 1:length(obj.DAP_table.Properties.VariableNames)
                name = obj.DAP_table.Properties.VariableNames{i};
                fprintf('%s\t',name)
                if rem(i,5) == 0
                    fprintf('\n')
                end
            end
        end

        function run_linear_model(obj)
            addpath(genpath('fdr_bh'))
            for i_task = 1:obj.n_task
                for i = 1:height(obj.T_ROI_position)
                    for j = 1:height(obj.T_ROI_position)
                        if i >= j
                            continue
                        end
    
                        name = [obj.task_flag{i_task}, '_', num2str(i), '_', num2str(j)];
                        set_formula(obj, name);
                        mdl = fitlm(obj.DAP_table,obj.formula);
                        
                        for i_var = 1:length(obj.independent_var)
                            ind_var = split(string(obj.independent_var{i_var}),":");
                            if length(ind_var) == 1
                                examine_var = mdl.CoefficientNames{contains(mdl.CoefficientNames, ind_var)};
                            else
                                examine_var = mdl.CoefficientNames{contains(mdl.CoefficientNames, ind_var(1))&...
                                    contains(mdl.CoefficientNames, ind_var(2))};
                            end
                            obj.model_table{i_var, [name,'_t']} = mdl.Coefficients{examine_var,3}; % t value
                            obj.model_table{i_var, [name,'_p']} = mdl.Coefficients{examine_var,4}; % p value
                            obj.model_table{i_var, [name,'_beta']} = mdl.Coefficients{examine_var,1}; % t value
                            obj.model_table{i_var, [name,'_n']} = mdl.NumObservations; % n value
                            obj.model_table{i_var, [name,'_SE']}   = mdl.Coefficients{examine_var,2}; % SE value
                        end
                    end
                end
            end
        end

        function add_multiple_comparision_control(obj)
            for i_task = 1:obj.n_task
                for i_var = 1:length(obj.independent_var)
                    model_table_names = obj.model_table.Properties.VariableNames;
                    p_index = find(endsWith(model_table_names,'_p') &...
                        startsWith(model_table_names,obj.task_flag{i_task}));
                    p_array = obj.model_table{i_var,p_index};
                    [~, ~, ~, p_array_BH]=fdr_bh(p_array,0.05,'pdep','no');
                    for i_name = 1:length(p_index)
                        obj.model_table{i_var, [model_table_names{p_index(i_name)},'_BH']} = p_array_BH(i_name); % p value
                    end
                end
            end
        end

        function set_formula(obj, name)
            obj.formula = [name,'~ 1'];
            for i = 1:length(obj.independent_var)
                obj.formula = [obj.formula, ' + ', obj.independent_var{i}];
            end
        end
        function plot_scatter(obj, task, ch1, ch2, var)
            figure
            x = obj.DAP_table.(var);
            name = [task, '_', num2str(ch1), '_', num2str(ch2)];
            y = obj.DAP_table.(name);
            plot(x, y, 'r.','MarkerSize',8)
            mdl = fitlm(x, y, 'linear');
            hold on;
            plot(x, predict(mdl, x), 'r', 'LineWidth', 2);
            i_var = find(strcmp(obj.independent_var,var));
            title(sprintf('roi%d - roi%d\nb=%.1f,t=%.1f,P=%.4f', ch1, ch2,...
                obj.model_table{i_var,[name,'_beta']}, ...
                obj.model_table{i_var,[name,'_t']},...
                obj.model_table{i_var,[name,'_p']}))
            ylabel('connectivity')
            xlabel(var)
            hold off;
        end
        function plot_connectivity(obj, threshold, clim_limit, FDR, var)
            i_var = contains(obj.independent_var, var);
            x = obj.T_ROI_position.x_mean;
            y = obj.T_ROI_position.y_mean;
            n_roi = height(obj.T_ROI_position);
            jet_color = jet(256);

            for i_task = 1:length(obj.task_flag)
                figure('name', [obj.task_flag{i_task}, '~', var])
                s = zeros((n_roi*n_roi-n_roi)/2,1);
                t = zeros((n_roi*n_roi-n_roi)/2,1);
                weights = zeros((n_roi*n_roi-n_roi)/2,1);
                k = 1;
                for i = 1:height(obj.T_ROI_position)
                    for j = 1:height(obj.T_ROI_position)
                        if i >= j
                            continue
                        end
                        s(k) = i; t(k) = j;
                        name = [obj.task_flag{i_task}, '_', num2str(i), '_', num2str(j)];
                        if FDR
                            weights(k) = obj.model_table{i_var, [name,'_p_BH']};
                        else
                            weights(k) = obj.model_table{i_var, [name,'_p']};
                        end
                        k = k + 1;
                    end
                end
                G = graph(s,t,weights);
                G = rmedge(G, find(abs(G.Edges{:,2})>threshold));
                G = rmedge(G, find(isnan(G.Edges{:,2})));
                
                if ~exist('clim_limit','var')
                    clim_limit = max(abs(G.Edges.Weight));
                end
                G.Edges.Weight(G.Edges.Weight>clim_limit) = clim_limit;
                G.Edges.Weight(G.Edges.Weight<-clim_limit) = -clim_limit;
                LWidths = G.Edges.Weight/clim_limit; % [0 1]
                color_array = jet_color(floor(LWidths*255)+1,:);
                
                plot(G,'XData',x,'YData',y,'edgecolor',color_array,'LineWidth',2)
                axis off
                colormap(jet_color) 
                colorbar
                % clim([0, clim_limit])
                set(gca, 'clim', [0 clim_limit]);
            end
        end

        function make_table(obj,varargin)
            i_var = find(strcmp(obj.independent_var, varargin{1}));
            T_all = table();
            for k = 1:obj.n_task
                ROI_1 = {};
                ROI_2 = {};
                fc_DAP_n = [];
                fc_DAP_beta = [];
                fc_DAP_SE = [];
                fc_DAP_t = [];
                fc_DAP_p = [];
                fc_DAP_p_BH = [];
                for i = 1:height(obj.T_ROI_position)
                    for j = 1:height(obj.T_ROI_position)
                        if i>=j
                            continue
                        end
%                         ROI_1(end+1) = obj.T_ROI_position{i,1};
%                         ROI_2(end+1) = obj.T_ROI_position{j,1};
                        ROI_1(end+1) = {i};
                        ROI_2(end+1) = {j};
                        
                        name = [obj.task_flag{k},'_',num2str(i),'_', num2str(j)];
                        
                        fc_DAP_n(end+1) = obj.model_table{i_var, [name,'_n']}; % n value
                        fc_DAP_beta(end+1)  = obj.model_table{i_var, [name,'_beta']}; % beta value
                        fc_DAP_SE(end+1)  = obj.model_table{i_var, [name,'_SE']}; % beta value
                        fc_DAP_t(end+1) = obj.model_table{i_var, [name,'_t']}; % t value
                        fc_DAP_p(end+1) = obj.model_table{i_var, [name,'_p']}; % p value
                        fc_DAP_p_BH(end+1) = obj.model_table{i_var, [name,'_p_BH']}; % p value
                    end
                end
                T = table(fc_DAP_n', fc_DAP_beta', fc_DAP_SE', fc_DAP_t',fc_DAP_p',fc_DAP_p_BH');
            
                T = renamevars(T,{'Var1','Var2','Var3','Var4','Var5','Var6'},...
                    {[obj.task_flag{k},'_n'],[obj.task_flag{k},'_beta'],...
                    [obj.task_flag{k},'_SE'],[obj.task_flag{k},'_t'],...
                    [obj.task_flag{k},'_p'],[obj.task_flag{k},'_p_BH']} );
                T_all = [T_all,T];
            end
            T_ROI = table(ROI_1',ROI_2');
            T_all = [T_ROI,T_all];
            T_all = renamevars(T_all,{'Var1','Var2'},{'ROI_1','ROI_2'} );
            if nargin > 1
                writetable(T_all, fullfile('..','table', ['sup_mat_',varargin{2},'.csv']))
            else
                writetable(T_all, fullfile('..','table', ['sup_mat_',varargin{1},'.csv']))
            end
        end

        function make_table_95CI(obj,varargin)
            % i haven't change anything here yet
            i_var = find(strcmp(obj.independent_var, varargin{1}));
            T_all = table();
            for k = 1:obj.n_task
                ROI_1 = {};
                ROI_2 = {};
                fc_DAP_n = [];
                fc_DAP_beta = [];
                fc_DAP_95CI = [];
                fc_DAP_beta_95CI = {};
                for i = 1:height(obj.T_ROI_position)
                    for j = 1:height(obj.T_ROI_position)
                        if i>=j
                            continue
                        end
                        ROI_1(end+1) = cellstr(sprintf('%s(%d)',char(obj.T_ROI_position{i,1}),i));
                        ROI_2(end+1) = cellstr(sprintf('%s(%d)',char(obj.T_ROI_position{j,1}),j));
                        
                        name = [obj.task_flag{k},'_',num2str(i),'_', num2str(j)];
                        
                        fc_DAP_n(end+1) = obj.model_table{i_var, [name,'_n']}; % n value
                        fc_DAP_beta(end+1)  = obj.model_table{i_var, [name,'_beta']}; % beta value
                        fc_DAP_95CI(end+1)   = 1.96*obj.model_table{i_var, [name,'_SE']};
                        fc_DAP_beta_95CI(end+1) = cellstr(sprintf('%.2f(%.2f,%.2f)',fc_DAP_beta(end), fc_DAP_beta(end)-fc_DAP_95CI(end),fc_DAP_beta(end)+fc_DAP_95CI(end)));
                    end
                end
                T = table(fc_DAP_n', fc_DAP_beta_95CI');
            
                T = renamevars(T,{'Var1','Var2'},...
                    {[obj.task_flag{k},'_n'],[obj.task_flag{k},'_beta_95CI']} );
                T_all = [T_all,T];
            end
            T_ROI = table(ROI_1',ROI_2');
            T_all = [T_ROI,T_all];
            T_all = renamevars(T_all,{'Var1','Var2'},{'ROI_1','ROI_2'} );
            if nargin > 1
                writetable(T_all, fullfile('..','table', ['sup_mat_',varargin{2},'.csv']))
            else
                writetable(T_all, fullfile('..','table', ['sup_mat_',varargin{1},'.csv']))
            end
        end

        function plot_range(obj,varargin)
            i_var = find(strcmp(obj.independent_var, varargin{1}));
            for k = 1:obj.n_task
                ROI_1 = {};
                ROI_2 = {};
                ROI_list = {};
                fc_DAP_beta = [];
                fc_DAP_SE = [];
                x = [];
                l = 1;
                ROI_list = {};
                for i = 1:height(obj.T_ROI_position)
                    for j = 1:height(obj.T_ROI_position)
                        if i>=j
                            continue
                        end
                        ROI_1(end+1) = obj.T_ROI_position{i,1};
                        ROI_2(end+1) = obj.T_ROI_position{j,1};
                        ROI_list{l} = [char(obj.T_ROI_position{i,1}),'-', char(obj.T_ROI_position{j,1})];
                        
                        name = [obj.task_flag{k},'_',num2str(i),'_', num2str(j)];
                        
                        fc_DAP_beta(end+1)  = obj.model_table{i_var, [name,'_beta']}; % beta value
                        fc_DAP_SE(end+1)  = obj.model_table{i_var, [name,'_SE']};
                        x(end+1) = l;
                        l = l + 1;
                    end
                end
                figure
                errorbar(fc_DAP_beta(1:53), x(1:53),fc_DAP_SE(1:53)*1.96,'horizontal','linestyle','none', 'marker','.', 'markersize',10,'color','blue');
                hold on
                plot([0 0],[0 120],'k-.')
                ylim([0 54])
                set(gca, 'YDir','reverse','Ycolor','blue')
                yticks(x(1:53))
                yticklabels(ROI_list(1:53))
                ax = gca;
                ax.XAxis.FontSize = 15;
                ax.YAxis.FontSize = 15;
                set(gcf, 'position', [10 10 835 1238])
                box off
                xlim([-0.6 0.6])
                xticks([-0.3 0.3])
                
                figure
                errorbar(fc_DAP_beta(54:end), x(54:end),fc_DAP_SE(54:end)*1.96,'horizontal','linestyle','none', 'marker','.', 'markersize',10,'color','red');
                hold on
                plot([0 0],[0 120],'k-.')
                ylim([53 106])
                set(gca, 'YDir','reverse','Ycolor','red')
                yticks(x(54:end))
                yticklabels(ROI_list(54:end))
                ax = gca;
                ax.XAxis.FontSize = 15;
                ax.YAxis.FontSize = 15;
                set(gca,'YAxisLocation','right')
                set(gcf, 'position', [10 10 973 1238])
                box off
                xlim([-0.6 0.6])
                xticks([-0.3 0.3])
            end
        end
        
        function add_behavior(obj)
            T = obj.DAP_table;
            T = convertvars(T,1,'string');
            T = convertvars(T,1,'double');
            % 
            % % wisconsin
            T_bx_old = readtable(fullfile('..','data', 'chm_all_bx_t.csv'));
            T_bx_old = T_bx_old(contains(T_bx_old{:,2}, 'card'),:);
            T_bx_old = T_bx_old(strcmp(T_bx_old{:,3}, 'match_perseverativeErrors'),:);
            T_bx_old = unstack(T_bx_old,'value','measure');
            T_bx_old{T_bx_old{:,"subj"} < 100000,1} = T_bx_old{T_bx_old{:,"subj"}<100000,1} * 100 + 4;
            T_bx_old(T_bx_old{:,"subj"} > 2196304,:) = [];
            
            T_bx_old(181,:) = []; % 286 subejcts
            T = outerjoin(T,T_bx_old, 'LeftKeys','id','RightKeys','subj','Type', 'left', 'MergeKeys', false);
            T.subj=[];
            T.task = [];
            % 
            % % word
            T_bx_old = readtable(fullfile('..','data', 'chm_all_bx_t.csv'));
            T_bx_old = T_bx_old(contains(T_bx_old{:,2}, 'word'),:);
            T_bx_old = T_bx_old(contains(T_bx_old{:,3}, 'meaning_acc')| contains(T_bx_old{:,3}, 'size_acc'),:);
            T_bx_old = unstack(T_bx_old,'value','measure');
            T_bx_old{T_bx_old{:,"subj"} < 100000,1} = T_bx_old{T_bx_old{:,"subj"}<100000,1} * 100 + 4;
            T_bx_old(T_bx_old{:,"subj"} > 2196304,:) = [];
            % u=unique(T_bx_old.subj);         % the unique values
            % [n,bin]=histc(T_bx_old.subj,u);  % count how many of each and where
            % ix1=find(n>1);  
            % % row 183 is repeated, 
            T_bx_old(183,:) = []; % 293 subejcts
            T_bx_old.meaning_size_acc = - T_bx_old.meaning_acc + T_bx_old.size_acc;%% because size acc should be higher than meaning acc
            T_bx_old(T_bx_old.meaning_acc == 0,:) = [];
            T = outerjoin(T,T_bx_old, 'LeftKeys','id','RightKeys','subj','Type', 'left', 'MergeKeys', false);
            T.subj=[];
            T.task = [];
            % 
            % % lmem
            T_bx_old = readtable('/Volumes/Projects/NIRS_Projects/ECHO/y18/derivatives_20230428_YY/all_bx_t.csv');
            T_bx_old = T_bx_old(strcmp(T_bx_old{:,2}, 'lmem'),:);
            % T_bx_old = T_bx_old(strcmp(T_bx_old{:,3}, 'acc'),:);
            T_bx_old = unstack(T_bx_old,'value','measure');
            T_bx_old{T_bx_old{:,"subj"} < 100000,1} = T_bx_old{T_bx_old{:,"subj"}<100000,1} * 100 + 4;
            T_bx_old(T_bx_old{:,"subj"} > 2196304,:) = [];
            % u=unique(T_bx_old.subj);         % the unique values
            % [n,bin]=histc(T_bx_old.subj,u);  % count how many of each and where
            % ix1=find(n>1);  
            % % row 183 is repeated, 
            T_bx_old(183,:) = []; % 293 subejcts
            T_bx_old(T_bx_old{:,"acc"} == 0,:) = [];
            
            T = outerjoin(T,T_bx_old, 'LeftKeys','id','RightKeys','subj','Type', 'left', 'MergeKeys', false);
            T.subj=[];
            T.task = [];
            obj.DAP_table = T;
        end
    end
end
