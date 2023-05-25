%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'CKmcmc';
M_.dynare_version = '5.3';
oo_.dynare_version = '5.3';
options_.dynare_version = '5.3';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'e1'};
M_.exo_names_tex(1) = {'e1'};
M_.exo_names_long(1) = {'e1'};
M_.exo_names(2) = {'e2'};
M_.exo_names_tex(2) = {'e2'};
M_.exo_names_long(2) = {'e2'};
M_.endo_names = cell(11,1);
M_.endo_names_tex = cell(11,1);
M_.endo_names_long = cell(11,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y'};
M_.endo_names_long(1) = {'y'};
M_.endo_names(2) = {'c'};
M_.endo_names_tex(2) = {'c'};
M_.endo_names_long(2) = {'c'};
M_.endo_names(3) = {'k'};
M_.endo_names_tex(3) = {'k'};
M_.endo_names_long(3) = {'k'};
M_.endo_names(4) = {'n'};
M_.endo_names_tex(4) = {'n'};
M_.endo_names_long(4) = {'n'};
M_.endo_names(5) = {'invest'};
M_.endo_names_tex(5) = {'invest'};
M_.endo_names_long(5) = {'invest'};
M_.endo_names(6) = {'z1'};
M_.endo_names_tex(6) = {'z1'};
M_.endo_names_long(6) = {'z1'};
M_.endo_names(7) = {'z2'};
M_.endo_names_tex(7) = {'z2'};
M_.endo_names_long(7) = {'z2'};
M_.endo_names(8) = {'MUC'};
M_.endo_names_tex(8) = {'MUC'};
M_.endo_names_long(8) = {'MUC'};
M_.endo_names(9) = {'MUL'};
M_.endo_names_tex(9) = {'MUL'};
M_.endo_names_long(9) = {'MUL'};
M_.endo_names(10) = {'r'};
M_.endo_names_tex(10) = {'r'};
M_.endo_names_long(10) = {'r'};
M_.endo_names(11) = {'w'};
M_.endo_names_tex(11) = {'w'};
M_.endo_names_long(11) = {'w'};
M_.endo_partitions = struct();
M_.param_names = cell(15,1);
M_.param_names_tex = cell(15,1);
M_.param_names_long = cell(15,1);
M_.param_names(1) = {'alppha'};
M_.param_names_tex(1) = {'alppha'};
M_.param_names_long(1) = {'alppha'};
M_.param_names(2) = {'betta'};
M_.param_names_tex(2) = {'betta'};
M_.param_names_long(2) = {'betta'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'delta'};
M_.param_names_long(3) = {'delta'};
M_.param_names(4) = {'gam'};
M_.param_names_tex(4) = {'gam'};
M_.param_names_long(4) = {'gam'};
M_.param_names(5) = {'nss'};
M_.param_names_tex(5) = {'nss'};
M_.param_names_long(5) = {'nss'};
M_.param_names(6) = {'rho1'};
M_.param_names_tex(6) = {'rho1'};
M_.param_names_long(6) = {'rho1'};
M_.param_names(7) = {'sigma1'};
M_.param_names_tex(7) = {'sigma1'};
M_.param_names_long(7) = {'sigma1'};
M_.param_names(8) = {'rho2'};
M_.param_names_tex(8) = {'rho2'};
M_.param_names_long(8) = {'rho2'};
M_.param_names(9) = {'sigma2'};
M_.param_names_tex(9) = {'sigma2'};
M_.param_names_long(9) = {'sigma2'};
M_.param_names(10) = {'psi'};
M_.param_names_tex(10) = {'psi'};
M_.param_names_long(10) = {'psi'};
M_.param_names(11) = {'c1'};
M_.param_names_tex(11) = {'c1'};
M_.param_names_long(11) = {'c1'};
M_.param_names(12) = {'iss'};
M_.param_names_tex(12) = {'iss'};
M_.param_names_long(12) = {'iss'};
M_.param_names(13) = {'yss'};
M_.param_names_tex(13) = {'yss'};
M_.param_names_long(13) = {'yss'};
M_.param_names(14) = {'kss'};
M_.param_names_tex(14) = {'kss'};
M_.param_names_long(14) = {'kss'};
M_.param_names(15) = {'css'};
M_.param_names_tex(15) = {'css'};
M_.param_names_long(15) = {'css'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 11;
M_.param_nbr = 15;
M_.orig_endo_nbr = 11;
M_.aux_vars = [];
options_.varobs = cell(2, 1);
options_.varobs(1)  = {'c'};
options_.varobs(2)  = {'n'};
options_.varobs_id = [ 2 4  ];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 11;
M_.eq_nbr = 11;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 4 0;
 0 5 0;
 1 6 0;
 0 7 0;
 0 8 0;
 2 9 0;
 3 10 0;
 0 11 15;
 0 12 0;
 0 13 16;
 0 14 0;]';
M_.nstatic = 6;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.dynamic_tmp_nbr = [4; 2; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'MUC' ;
  2 , 'name' , 'MUL' ;
  3 , 'name' , 'r' ;
  4 , 'name' , 'w' ;
  5 , 'name' , '5' ;
  6 , 'name' , '6' ;
  7 , 'name' , 'z1' ;
  8 , 'name' , 'z2' ;
  9 , 'name' , 'y' ;
  10 , 'name' , 'invest' ;
  11 , 'name' , 'k' ;
};
M_.mapping.y.eqidx = [9 10 ];
M_.mapping.c.eqidx = [1 10 ];
M_.mapping.k.eqidx = [3 4 9 11 ];
M_.mapping.n.eqidx = [3 4 9 ];
M_.mapping.invest.eqidx = [10 11 ];
M_.mapping.z1.eqidx = [3 4 7 9 ];
M_.mapping.z2.eqidx = [2 8 ];
M_.mapping.MUC.eqidx = [1 5 6 ];
M_.mapping.MUL.eqidx = [2 6 ];
M_.mapping.r.eqidx = [3 5 ];
M_.mapping.w.eqidx = [4 6 ];
M_.mapping.e1.eqidx = [7 ];
M_.mapping.e2.eqidx = [8 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 6 7 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(11, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(15, 1);
M_.endo_trends = struct('deflator', cell(11, 1), 'log_deflator', cell(11, 1), 'growth_factor', cell(11, 1), 'log_growth_factor', cell(11, 1));
M_.NNZDerivatives = [34; -1; -1; ];
M_.static_tmp_nbr = [4; 2; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
close all;
set_param_value('alppha',0.33)
set_param_value('betta', 0.99)
set_param_value('delta',0.025)
set_param_value('gam', 2.0)
set_param_value('rho1', 0.9)
set_param_value('sigma1', 0.02)
set_param_value('rho2', 0.7)
set_param_value('sigma2', 0.01)
set_param_value('nss', 1/3)
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
estim_params_.var_exo = zeros(0, 10);
estim_params_.var_endo = zeros(0, 10);
estim_params_.corrx = zeros(0, 11);
estim_params_.corrn = zeros(0, 11);
estim_params_.param_vals = zeros(0, 10);
estim_params_.param_vals = [estim_params_.param_vals; 2, NaN, (-Inf), Inf, 5, NaN, NaN, 0.95, 0.9999, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, NaN, (-Inf), Inf, 5, NaN, NaN, 0.0, 5.0, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 6, NaN, (-Inf), Inf, 5, NaN, NaN, 0.0, 1.0, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 7, NaN, (-Inf), Inf, 5, NaN, NaN, 0.0, 0.1, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 8, NaN, (-Inf), Inf, 5, NaN, NaN, 0.0, 1.0, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 9, NaN, (-Inf), Inf, 5, NaN, NaN, 0.0, 0.1, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 5, NaN, (-Inf), Inf, 5, NaN, NaN, 0.25, 0.375, NaN ];
options_.irf = 40;
options_.mh_init_scale = 5;
options_.mh_jscale = 0.8;
options_.mh_nblck = 1;
options_.mh_replic = 10000;
options_.order = 1;
options_.datafile = 'dsgedata';
options_.nobs = 160;
var_list_ = {};
oo_recursive_=dynare_estimation(var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'CKmcmc_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'CKmcmc_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'CKmcmc_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'CKmcmc_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'CKmcmc_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'CKmcmc_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'CKmcmc_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
