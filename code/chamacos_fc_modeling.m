close all
clc
clear all
%% add folder fdr_bh to the path
addpath('fdr_bh')
%% load connectivity data
load(fullfile('..','data','fc_card_roi.mat'), 'fc_card');
load(fullfile('..','data','fc_stern_roi.mat'), 'fc_stern');
load(fullfile('..','data','fc_ppt_roi.mat'), 'fc_ppt');

%% init
fc = Connectivity();
fc.set_DAP_table;
fc.set_T_ROI_position;
%% add task
fc.add_task(fc_card.match-fc_card.control, 'card');
fc.add_task(fc_stern.lmem, 'lmem');
fc.add_task(fc_ppt.word - fc_ppt.size, 'ppt');

%% add variables
warning off
independent_var = {'ldap_tot_i_sg2_pg', 'age_fnirs', ...
    'sex', 'momdl_age2','educcat_mom', 'z2home_6','usyr3', ...
    'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model();
fc.add_behavior;DAP_conn_behav = fc.DAP_table;
writetable(DAP_conn_behav,fullfile('..','data','DAP_conn_behav.csv'))
fc.add_multiple_comparision_control;
independent_var = {'ldap_tot_i_sg2_pg', 'age_fnirs', ...
    'sex', 'momdl_age2','educcat_mom', 'z2home_6','usyr3', ...
    'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model;
fc.add_multiple_comparision_control;
save(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
writetable(fc.model_table, fullfile('..','processed_data', 'model_table_preg.csv'))

independent_var = {'ldap_tot_c_auc', 'age_fnirs', ...
    'sex', 'momdl_age2','educcat_mom', 'z2home_6','usyr3', ...
    'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model;
fc.add_multiple_comparision_control;
save(fullfile('..','processed_data', 'fc_child.mat'), 'fc')
writetable(fc.model_table, fullfile('..','processed_data', 'model_table_child.csv'))

%% all participants, DAP plot
load(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
fc.plot_connectivity(0.05, 0.05, true, 'ldap_tot_i_sg2_pg')
load(fullfile('..','processed_data', 'fc_child.mat'), 'fc')
fc.plot_connectivity(0.05, 0.05, true, 'ldap_tot_c_auc')

%% all participants, DAP, by sex

load(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
index = find(fc.DAP_table.sex == 'Boy');
fc.stratify(index);
independent_var = {'ldap_tot_i_sg2_pg', 'age_fnirs', 'momdl_age2',...
    'educcat_mom', 'z2home_6','usyr3', 'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model();
fc.add_multiple_comparision_control;
save(fullfile('..','processed_data', 'fc_preg_boy.mat'), 'fc')

load(fullfile('..','processed_data', 'fc_child.mat'), 'fc')
index = find(fc.DAP_table.sex == 'Boy');
fc.stratify(index);
independent_var = {'ldap_tot_c_auc', 'age_fnirs', 'momdl_age2',...
    'educcat_mom', 'z2home_6','usyr3', 'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model();
fc.add_multiple_comparision_control;
save(fullfile('..','processed_data', 'fc_child_boy.mat'), 'fc')


load(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
index = find(fc.DAP_table.sex == 'Girl');
fc.stratify(index);
independent_var = {'ldap_tot_i_sg2_pg', 'age_fnirs', 'momdl_age2',...
    'educcat_mom', 'z2home_6','usyr3', 'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model();
fc.add_multiple_comparision_control;
save(fullfile('..','processed_data', 'fc_preg_girl.mat'), 'fc')

load(fullfile('..','processed_data', 'fc_child.mat'), 'fc')
index = find(fc.DAP_table.sex == 'Girl');
fc.stratify(index);
independent_var = {'ldap_tot_c_auc', 'age_fnirs', 'momdl_age2',...
    'educcat_mom', 'z2home_6','usyr3', 'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model();
fc.add_multiple_comparision_control;
save(fullfile('..','processed_data', 'fc_child_girl.mat'), 'fc')


load(fullfile('..','processed_data', 'fc__preg_boy.mat'), 'fc')
fc.plot_connectivity(0.05, 0.05, true, 'ldap_tot_i_sg2_pg')
load(fullfile('..','processed_data', 'fc_child_boy.mat'), 'fc')
fc.plot_connectivity(0.05, 0.05, true, 'ldap_tot_c_auc')


load(fullfile('..','processed_data', 'fc_preg_girl.mat'), 'fc')
fc.plot_connectivity(0.05, 0.05, true, 'ldap_tot_i_sg2_pg')
load(fullfile('..','processed_data', 'fc_child_girl.mat'), 'fc')
fc.plot_connectivity(0.05, 0.05, true, 'ldap_tot_c_auc')

%% all participants, DAP, affected by cov
load(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
fc.plot_connectivity(0.05, 0.05, true, 'sex')
load(fullfile('..','processed_data', 'fc_child.mat'), 'fc')
fc.plot_connectivity(0.05, 0.05, true, 'sex')

%% make tables for sup material
load(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
fc.make_table('ldap_tot_i_sg2_pg','sep_preg_95CI')
load(fullfile('..','processed_data', 'fc_child.mat'), 'fc')
fc.make_table('ldap_tot_c_auc','sep_child_95CI')

load(fullfile('..','processed_data', 'fc__preg_boy.mat'), 'fc')
fc.make_table('ldap_tot_i_sg2_pg','sep_preg_boy')
load(fullfile('..','processed_data', 'fc_child_boy.mat'), 'fc')
fc.make_table('ldap_tot_c_auc','sep_child_boy')
load(fullfile('..','processed_data', 'fc_preg_girl.mat'), 'fc')
fc.make_table('ldap_tot_i_sg2_pg','sep_preg_girl')
load(fullfile('..','processed_data', 'fc_child_girl.mat'), 'fc')
fc.make_table('ldap_tot_c_auc','sep_child_girl')

load(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
fc.make_table_95CI('ldap_tot_i_sg2_pg','sep_preg_95CI')
load(fullfile('..','processed_data', 'fc_child.mat'), 'fc')
fc.make_table_95CI('ldap_tot_c_auc','sep_child_95CI')

load(fullfile('..','processed_data', 'fc__preg_boy.mat'), 'fc')
fc.make_table_95CI('ldap_tot_i_sg2_pg','sep_preg_boy_95CI')
load(fullfile('..','processed_data', 'fc_child_boy.mat'), 'fc')
fc.make_table_95CI('ldap_tot_c_auc','sep_child_boy_95CI')
load(fullfile('..','processed_data', 'fc_preg_girl.mat'), 'fc')
fc.make_table_95CI('ldap_tot_i_sg2_pg','sep_preg_girl_95CI')
load(fullfile('..','processed_data', 'fc_child_girl.mat'), 'fc')
fc.make_table_95CI('ldap_tot_c_auc','sep_child_girl_95CI')

%% plot range
% load(fullfile('..','processed_data', 'fc.mat'), 'fc')
% fc.plot_range('ldap_tot_i_sg2_pg')
% fc.plot_range('ldap_tot_c_auc')

load(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
fc.plot_range('ldap_tot_i_sg2_pg')
load(fullfile('..','processed_data', 'fc_child.mat'), 'fc')
fc.plot_range('ldap_tot_c_auc')
%% include the interaction term sex and DAPs
load(fullfile('..','processed_data', 'fc_preg.mat'), 'fc')
independent_var = {'ldap_tot_i_sg2_pg', 'sex:ldap_tot_i_sg2_pg', 'age_fnirs', ...
    'sex', 'momdl_age2','educcat_mom', 'z2home_6','usyr3', ...
    'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model;
fc.add_multiple_comparision_control;
fc.make_table('sex:ldap_tot_i_sg2_pg','sep_preg_sex')
fc.make_table('ldap_tot_i_sg2_pg','sep_preg_sex_dap')

independent_var = {'ldap_tot_c_auc', 'sex:ldap_tot_c_auc', 'age_fnirs', ...
    'sex', 'momdl_age2','educcat_mom', 'z2home_6','usyr3', ...
    'povcat2_m18y_i','sppstd_i'};
fc.independent_var = independent_var;
fc.run_linear_model;
fc.add_multiple_comparision_control;
fc.make_table('sex:ldap_tot_c_auc','sep_child_sex')
fc.make_table('ldap_tot_c_auc','sep_child_sex_dap')