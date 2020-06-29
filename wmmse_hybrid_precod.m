function [obj_fun_reg_wmmse,vec_coupling,sum_rate_PDD,sum_rate,SpEff_final,beam_var] = wmmse_hybrid_precod(par,H_UL_effec,H_DL_effec,H_SI) %#ok<*INUSL,*INUSD>
%wmmse_hybrid_precod minimizes the WMMSE in order to maximize the sum rate.
%The regularization from Shi2017 is used as a means to overcome the analog
%modulus constraints
% INPUT
%  par             -- Struct with predefined parameters
%  H_UL_effec      -- Matrix gain between users and BS in UL
%  H_DL_effec      -- Matrix gain between users and BS in DL
%  Hmm_effec       -- Matrix gain between each user
%  H_SI            -- Matrix of SI gains
% OUTPUT
% QD_BF            -- Downlink beamforming matrices in a cell
% sta              -- Uplink beamforming values in a vector

%% Assign initial values for all major variables

% Baseband combiners - may not need to be initialized
% Get UL random vector
vec_qu_BB = sqrt(1/2)*complex(randn(par.antBS_RF,par.lambdaul),randn(par.antBS_RF,par.lambdaul));
% Get DL random scalar
scal_vd_BB = sqrt(1/2)*complex(randn(par.antUE,par.lambdaul),randn(par.antUE,par.lambdaul));

% Baseband precoders
% Get UL random scalar
scal_wu_BB = sqrt(1/2)*complex(randn(par.antUE,par.lambdaul),randn(par.antUE,par.lambdaul));
% Get DL random vector
vec_fbs_BB = sqrt(1/2)*complex(randn(par.antBS_RF,par.lambdaul),randn(par.antBS_RF,par.lambdaul));

% Auxiliary variables
% UL
scal_zu_aux = sqrt(1/2)*complex(randn(par.antUE,par.lambdaul),randn(par.antUE,par.lambdaul));
% DL
vec_zd_aux = sqrt(1/2)*complex(randn(par.antBS_Tx,par.lambdadl),randn(par.antBS_Tx,par.lambdadl));

% Necessary variables for the analog combiners/precoders depend on the
% phase shifter resolution
if strcmp(par.JointBF,'LOWRES_PHASE') % If finite resolution
    % Set of values for the analog combiners and precoders
    % Quantized set of UL users
    quant_set_ul = exp((2*pi*1j/2^(par.bit_res))*(1:2^(par.bit_res)))/sqrt(par.antUE);
    % Quantized set of DL users
    quant_set_dl = exp((2*pi*1j/2^(par.bit_res))*(1:2^(par.bit_res)));%/sqrt(par.antBS_Tx);
    
    % Analog combiners
    % UL matrix
    mat_Qrf_ul = zeros(par.antBS_Rx,par.antBS_RF);
    rand_ul_matrix = randi(2^(par.bit_res),par.antBS_Rx,par.antBS_RF);
    for idxRow = 1:size(mat_Qrf_ul,1)
        for idxCol = 1:size(mat_Qrf_ul,2)
            mat_Qrf_ul(idxRow,idxCol) = quant_set_ul(rand_ul_matrix(idxRow,idxCol));
        end
    end
    % DL scalar
    scal_vrf_dl = quant_set_dl(randi(2^(par.bit_res)));
    
    % Analog precoders
    % UL scalar
    scal_wrf_ul = quant_set_ul(randi(2^(par.bit_res)));
    % DL matrix
    mat_Frf_dl = zeros(par.antBS_Tx,par.antBS_RF);
    rand_dl_matrix = randi(2^(par.bit_res),par.antBS_Tx,par.antBS_RF);
    for idxRow = 1:size(mat_Frf_dl,1)
        for idxCol = 1:size(mat_Frf_dl,2)
            mat_Frf_dl(idxRow,idxCol) = quant_set_dl(rand_dl_matrix(idxRow,idxCol));
        end
    end
    clear rand_ul_matrix rand_dl_matrix;
    % Create the set to update the blocks in the BCD - define max and min
    idx_rand_min = 1;
    idx_rand_max = 6;%12
    
else % If infinite resolution
    % Set of values for the analog combiners and precoders
    % Quantized set of UL and DL users is the complex set, so we do not
    % define it
    
    % Analog combiners
    % UL matrix
    mat_Qrf_ul = complex(randn(par.antBS_Rx,par.antBS_RF),randn(par.antBS_Rx,par.antBS_RF));
    for idxRow = 1:size(mat_Qrf_ul,1)
        for idxCol = 1:size(mat_Qrf_ul,2)
            % Unitary components
            mat_Qrf_ul(idxRow,idxCol) = mat_Qrf_ul(idxRow,idxCol)/(abs(mat_Qrf_ul(idxRow,idxCol))*...
                sqrt(par.antUE));
        end
    end
    % DL scalar
    scal_vrf_dl = complex(randn,randn);
    scal_vrf_dl = scal_vrf_dl/(abs(scal_vrf_dl));%(abs(scal_vrf_dl)*sqrt(par.antBS_Tx));
    
    % Analog precoders
    % UL scalar
    scal_wrf_ul =  complex(randn,randn);
    scal_wrf_ul =  scal_wrf_ul/(abs(scal_wrf_ul)*sqrt(par.antUE));
    % DL matrix
    mat_Frf_dl = complex(randn(par.antBS_Tx,par.antBS_RF),randn(par.antBS_Tx,par.antBS_RF));
    for idxRow = 1:size(mat_Frf_dl,1)
        for idxCol = 1:size(mat_Frf_dl,2)
            mat_Frf_dl(idxRow,idxCol) = mat_Frf_dl(idxRow,idxCol)/...
                (abs(mat_Frf_dl(idxRow,idxCol)));%(abs(mat_Frf_dl(idxRow,idxCol))*sqrt(par.antBS_Tx));
        end
    end
    % Create the set to update the blocks in the BCD - define max and min
    idx_rand_min = -2;%-4;
    idx_rand_max = 4;%8;
end

% Weights for WMMSE
% Unitary weights for UL and DL
rho_weight_DL_ite = ones(par.lambdadl,1);
rho_weight_UL_ite = ones(par.lambdaul,1);

% Channel in the DL should have a transposed dimension, but not Hermitian
H_DL_effec = transp(H_DL_effec);

%% Definitions for PDD
% Constant values for stopping criteria
control_const = 0.8;
eta_0 = 1e-3;
epsilon_0 = 1e-3; %#ok<*NASGU>
episilon_zero = 1e-6;
max_iter_BCD = 30;
max_iter_Penalty = 50;
% Penalty parameters
delta_unique = zeros(max_iter_BCD,1);
delta_unique(1) = 100/par.antBS_Tx;%2/64;%par.antBS_Tx;% proposed by Shi2018 is 100/par.antBS;
% Initial value for Lagrangian multipliers in the UL and DL
scal_lambda_ul = zeros(max_iter_BCD,1);
vec_lambda_dl = zeros(par.antBS_Tx,max_iter_BCD);

% Vectors for mse and rate
% MSE
mse_ul = zeros(max_iter_BCD,1);
mse_dl = zeros(max_iter_BCD,1);
% Rate for BCD
rate_ul = zeros(max_iter_BCD,1);
rate_dl = zeros(max_iter_BCD,1);
sum_rate = zeros(max_iter_BCD,1);
% Rate for PDD
sum_rate_PDD = zeros(max_iter_Penalty,1);

% Vectors for convergence related functions
% Objective function with regularization
obj_fun_reg_wmmse = zeros(max_iter_BCD,1);
% Objective function for BCD
obj_fun_wmmse = zeros(max_iter_BCD,1);
% Objective function for PDD
obj_fun_wmmse_PDD = zeros(max_iter_Penalty,1);
% Augmented Lagrangian function
aug_Lagrange = zeros(max_iter_BCD,1);
% Relative error
relative_error = zeros(max_iter_BCD,1);
% Variables to keep track of changes within blocks
wmse_eval_iter = zeros(10,length(idx_rand_min:idx_rand_max));
rate_eval_iter = zeros(10,length(idx_rand_min:idx_rand_max));
iter_ord_blk = zeros(10,length(idx_rand_min:idx_rand_max));


% Coupling constraints may define the stopping criteria for UL and DL
coupling_h_ul_norm_two = zeros(max_iter_Penalty,1);
coupling_h_ul_norm_two(1) = norm(scal_zu_aux - scal_wrf_ul*scal_wu_BB);
coupling_h_ul_norm_inf = zeros(max_iter_Penalty,1);
coupling_h_dl_norm_two = zeros(max_iter_Penalty,1);
coupling_h_dl_norm_two(1) = norm((vec_zd_aux - mat_Frf_dl*vec_fbs_BB));
coupling_h_dl_norm_inf = zeros(max_iter_Penalty,1);

% Iteration of the BCD method
iter_BCD = 1;
iter_Penalty = 1;

% Continue until termination criteria is met
while (iter_Penalty <= max_iter_Penalty)
    %% Solve the sBSUM algorithm from [Shi2017]
    % First loop to solve the problem for a fixed Lagrangian multipliers
    % and penalties
    sBSUM_hybrid_iter;%sBSUM_hybrid_iter;_gsplit
    
    % Update the objective function and rate for each PDD iteration
    obj_fun_wmmse_PDD(iter_Penalty) = obj_fun_wmmse(iter_BCD-1);
    sum_rate_PDD(iter_Penalty) = sum_rate(iter_BCD-1);
    
    % Update coupling constraints with infinite norm
    % UL
    coupling_h_ul_norm_inf(iter_Penalty) = max(abs(scal_zu_aux - scal_wrf_ul*scal_wu_BB) );
    
    % DL
    coupling_h_dl_norm_inf(iter_Penalty) = max(abs(vec_zd_aux - mat_Frf_dl*vec_fbs_BB));
    
    %% Update Lagrangian multipliers
    % If regularization norm is lower than threshold, update Lagrange
    % multiplier for both
    if max([coupling_h_ul_norm_inf(iter_Penalty);coupling_h_dl_norm_inf(iter_Penalty)]) < eta_0
        % Update only lambda
        % UL
        scal_lambda_ul(iter_Penalty+1) = scal_lambda_ul(iter_Penalty) +...
            real(scal_zu_aux - scal_wrf_ul*scal_wu_BB)/delta_unique(iter_Penalty);
        % DL
        vec_lambda_dl(:,iter_Penalty+1) = vec_lambda_dl(:,iter_Penalty) +...
            real((vec_zd_aux - mat_Frf_dl*vec_fbs_BB)/delta_unique(iter_Penalty));
        
        % Penalty remains the same
        delta_unique(iter_Penalty + 1) = delta_unique(iter_Penalty);
        % Also update eta_0
%         eta_0 = control_const*min([eta_0 max([coupling_h_ul_norm_inf(iter_Penalty);coupling_h_dl_norm_inf(iter_Penalty)])]);
        eta_0 = 0.9*min([coupling_h_ul_norm_inf(iter_Penalty);coupling_h_dl_norm_inf(iter_Penalty)]);
    else
        % Lambda remains the same
        % UL
        scal_lambda_ul(iter_Penalty+1) = scal_lambda_ul(iter_Penalty);
        % DL
        vec_lambda_dl(:,iter_Penalty+1) = vec_lambda_dl(:,iter_Penalty);
        
        % Update the penalty
        delta_unique(iter_Penalty + 1) = control_const*delta_unique(iter_Penalty);
    end
    
    %% Evaluate the termination criteria
    if (max( [coupling_h_ul_norm_inf(iter_Penalty) coupling_h_dl_norm_inf(iter_Penalty)] ) < episilon_zero) && iter_Penalty >= 5
        break;
    end
%     if mod(iter_Penalty,5) == 0
%         disp('stop');
%     end
    % Update the value of the convergence parameter epsilon_0
    epsilon_0 = control_const*epsilon_0;
    % Increase Penalty iteration
    iter_Penalty = iter_Penalty + 1;
    % Update the maximum number of iterations if it converged before
    max_iter_BCD = iter_BCD - 1 + max_iter_Penalty;
end

%% Return variables
% Beamforming variables and gains to return
% UL
beam_var.UL.vec_qu_BB = vec_qu_BB;
beam_var.UL.mat_Qrf_ul = mat_Qrf_ul;
beam_var.UL.scal_wrf_ul = scal_wrf_ul;
beam_var.UL.scal_wu_BB = scal_wu_BB;
beam_var.UL.beam_gain_ul_tx = (abs((vec_qu_BB'*mat_Qrf_ul')*H_UL_effec*(scal_wrf_ul*scal_wu_BB) )^2);
beam_var.UL.beam_gain_ul_rx = real((vec_qu_BB'*mat_Qrf_ul')*mat_Psi_ul_upd*(mat_Qrf_ul*vec_qu_BB));
beam_var.UL.H_UL_effec = H_UL_effec;
% DL
beam_var.DL.scal_vd_BB = scal_vd_BB;
beam_var.DL.scal_vrf_dl = scal_vrf_dl;
beam_var.DL.mat_Frf_dl = mat_Frf_dl;
beam_var.DL.vec_fbs_BB = vec_fbs_BB;
beam_var.DL.beam_gain_dl_tx = ( abs((scal_vd_BB'*scal_vrf_dl')*H_DL_effec'*(mat_Frf_dl*vec_fbs_BB) )^2 );
beam_var.DL.beam_gain_dl_rx = real( (abs(scal_vrf_dl*scal_vd_BB)^2)*scal_Psi_dl );
beam_var.DL.H_DL_effec = H_DL_effec;

% Removing exceding dimensions before convergence
obj_fun_wmmse_PDD(iter_Penalty:end) = [];
vec_coupling = [coupling_h_ul_norm_inf coupling_h_dl_norm_inf];
sum_rate(iter_BCD:end) = [];
SpEff_final = [rate_ul(iter_BCD-1);rate_dl(iter_BCD-1)];
end