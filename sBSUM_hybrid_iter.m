%% File to include the sBSUM solution
% Corresponds to the first loop to solve the problem for fixed Lagrangian
% multipliers

% Iterate until convergence
while (iter_BCD <= max_iter_BCD)
    % Generate the sequence of random updates - 6 blocks
%     fst_rand = randi([idx_rand_min idx_rand_max],1);
%     while fst_rand == 0
%         fst_rand = randi([idx_rand_min idx_rand_max],1);
%     end
%     
%     % Order of updates
%     ord_upd = horzcat(fst_rand,idx_rand_min:fst_rand-1, fst_rand+1:idx_rand_max);%idx_rand_min:idx_rand_max;
    ord_upd = datasample(idx_rand_min:idx_rand_max,length(idx_rand_min:idx_rand_max),'Replace',false);
    ord_upd(ord_upd == 0) = [];
    
%     % Index to keep track of appearance in the vector
%     idxBlk_upd = 1;
%     iter_ord_blk(iter_BCD,:) = ord_upd;
    
    %% Update the interference + noise covariance
    % Evaluate the covariance matrix and scalar of the interference plus
    % noise
    mat_Psi_ul = H_SI*(vec_zd_aux*vec_zd_aux')*H_SI' + par.noise*eye(par.antBS_Rx);
    scal_Psi_dl = par.noise; % No user-to-user interference
    
    %% Loop in the vector to update this vector
    for idxLoop = ord_upd
        switch idxLoop
            
            case 1 % Baseband combiners
                %% Solving the baseband combiners block
                
                % Update the UL baseband combiner
                vec_qu_BB = pinv(abs(scal_zu_aux)^2*mat_Qrf_ul'*(H_UL_effec*H_UL_effec')*mat_Qrf_ul + ...
                    mat_Qrf_ul'*mat_Psi_ul*mat_Qrf_ul)*(mat_Qrf_ul'*H_UL_effec*scal_zu_aux);
                % Update the DL baseband combiner
                scal_vd_BB = (scal_vrf_dl'*H_DL_effec'*vec_zd_aux)/(H_DL_effec'*(vec_zd_aux*vec_zd_aux')*H_DL_effec ...
                    + scal_Psi_dl);
                
%                 % Update the WMMSE and sum rate for the block
%                 wmse_short_eval;
%                 % Saving the iteration for this block
%                 wmse_eval_iter(iter_BCD,idxBlk_upd) = sum_wmse_total_blk;
%                 rate_eval_iter(iter_BCD,idxBlk_upd) = sum_rate_blk;
%                 % Increase block count
%                 idxBlk_upd = idxBlk_upd + 1;
                
            case 2 % Baseband precoders
                %% Baseband precoders block
                % UL baseband precoder
                scal_wu_BB = (scal_zu_aux + delta_unique(iter_Penalty)*scal_lambda_ul(iter_Penalty))/scal_wrf_ul;
                % Get DL random vector
                vec_fbs_BB = pinv(mat_Frf_dl)*(vec_zd_aux + delta_unique(iter_Penalty)*vec_lambda_dl(:,iter_Penalty));
                
%                 % Update the WMMSE and sum rate for the block
%                 wmse_short_eval;
%                 % Saving the iteration for this block
%                 wmse_eval_iter(iter_BCD,idxBlk_upd) = sum_wmse_total_blk;
%                 rate_eval_iter(iter_BCD,idxBlk_upd) = sum_rate_blk;
%                 % Increase block count
%                 idxBlk_upd = idxBlk_upd + 1;
                
            case 3 % Auxiliary variables
                %% Auxiliary variables
                % UL auxiliary variable - scalar
                % Defining M_u and l_u - no user to user interference for now
                scal_Mu = rho_weight_UL_ite*H_UL_effec'*mat_Qrf_ul*(vec_qu_BB*vec_qu_BB')*mat_Qrf_ul'*H_UL_effec + ...
                    1/(2*delta_unique(iter_Penalty));
                scal_lu = rho_weight_UL_ite*H_UL_effec'*mat_Qrf_ul*vec_qu_BB + ...
                    (scal_wrf_ul*scal_wu_BB - delta_unique(iter_Penalty)*scal_lambda_ul(iter_Penalty))/(2*delta_unique(iter_Penalty));
                
                % Solve the problem initially with CVX
                cvx_begin quiet
                % Defines variables
                variable scal_Zu(1) complex;
                % Defining expression
                expression obj_fun;
                obj_fun = 0;
                obj_fun = obj_fun + real(quad_form(scal_Zu,scal_Mu)) -2*real(scal_Zu'*scal_lu);
                % Constraint related to the thrust region
                scal_Zu*scal_Zu' <= par.pmaxUL; %#ok<*VUNUS>
                % Defining the optimization problem
                minimize( obj_fun )
                cvx_end
                % Update the UL auxiliary variable
                scal_zu_aux = scal_Zu;
                
%                 auxiliar_zu_asympt;
                
                if ~par.asympt % no asymptotic analysis in the DL
                    % DL auxiliary variable - vector
                    % Defining M_d and l_d
                    mat_Md = rho_weight_UL_ite*H_SI'*mat_Qrf_ul*(vec_qu_BB*vec_qu_BB')*mat_Qrf_ul'*H_SI + ...
                        rho_weight_DL_ite*abs(scal_vrf_dl*scal_vd_BB)^2*(H_DL_effec*H_DL_effec') + ...
                        (1/(2*delta_unique(iter_Penalty)))*eye(par.antBS_Tx);
                    vec_ld = rho_weight_DL_ite*(scal_vrf_dl*scal_vd_BB)*H_DL_effec + (1/(2*delta_unique(iter_Penalty)))*...
                        (mat_Frf_dl*vec_fbs_BB - delta_unique(iter_Penalty)*vec_lambda_dl(:,iter_Penalty) );
                    
                    % Solve the problem initially with CVX
                    cvx_begin quiet
                    % Defines variables
                    variable vec_Zd(par.antBS_Tx,1) complex;
                    dual variable phi_lagr_DL;
                    % Defining expression
                    expression obj_fun;
                    obj_fun = 0;
                    try
                        obj_fun = obj_fun + real(quad_form(vec_Zd,mat_Md)) -2*real(vec_Zd'*vec_ld);
                    catch
                        mat_Md = nearestSPD(mat_Md);
                        obj_fun = obj_fun + real(quad_form(vec_Zd,mat_Md)) -2*real(vec_Zd'*vec_ld);
                        str_save = strcat('Iter_BCD',num2str(iter_BCD));
                        % Save that this iteration had a psd problem
                        save(strcat(par.sta_save_name,'_',num2str(par.seed),'_psd'),'str_save');
                    end
                    % Constraint related to the thrust region
                    phi_lagr_DL : real(quad_form(vec_Zd,eye(par.antBS_Tx))) <= par.pmaxDL;
                    % Defining the optimization problem
                    minimize( obj_fun )
                    cvx_end
                    % Update the auxiliary variable
                    vec_zd_aux = vec_Zd;
                
                else % with asymptotic analysis in the DL
                    auxiliar_zd_asympt;
                end
                
                % Evaluate the covariance matrix and scalar of the interference plus
                % noise - Updated after this iteration because z_d is
                % present here
                mat_Psi_ul = H_SI*(vec_zd_aux*vec_zd_aux')*H_SI' + par.noise*eye(par.antBS_Rx);
                scal_Psi_dl = par.noise; % No user-to-user interference
                
%                 % Update the WMMSE and sum rate for the block
%                 wmse_short_eval;
%                 % Saving the iteration for this block
%                 wmse_eval_iter(iter_BCD,idxBlk_upd) = sum_wmse_total_blk;
%                 rate_eval_iter(iter_BCD,idxBlk_upd) = sum_rate_blk;
%                 % Increase block count
%                 idxBlk_upd = idxBlk_upd + 1;
                
            case 4 % WMMSE weights
                %% Update the WMMSE weights
                % Evaluate the UL MSE
                mse_ul(iter_BCD) = real(abs(1 - (vec_qu_BB'*mat_Qrf_ul')*H_UL_effec*(scal_zu_aux) )^2 + ...
                    (vec_qu_BB'*mat_Qrf_ul')*mat_Psi_ul*(mat_Qrf_ul*vec_qu_BB)); %#ok<*SAGROW>
                % Update the UL weights
                rho_weight_UL_ite = 1/mse_ul(iter_BCD);
                
                % Evaluate the DL MSE
                mse_dl(iter_BCD) = real(abs(1 - (scal_vd_BB'*scal_vrf_dl')*H_DL_effec'*vec_zd_aux )^2 + ...
                    (abs(scal_vrf_dl*scal_vd_BB)^2)*scal_Psi_dl);
                % Update the DL weights
                rho_weight_DL_ite = 1/mse_dl(iter_BCD);
                
%                 % Update the WMMSE and sum rate for the block
%                 wmse_short_eval;
%                 % Saving the iteration for this block
%                 wmse_eval_iter(iter_BCD,idxBlk_upd) = sum_wmse_total_blk;
%                 rate_eval_iter(iter_BCD,idxBlk_upd) = sum_rate_blk;
%                 % Increase block count
%                 idxBlk_upd = idxBlk_upd + 1;
                
            case 5 % Analog combiners for finite resolution
                %% Analog Combiners Block
                % UL analog combiner - matrix
                % Define matrices J_Q, K_Q, and N_Q
                mat_J_Q = rho_weight_UL_ite*(abs(scal_zu_aux)^2*(H_UL_effec*H_UL_effec') + mat_Psi_ul );
                mat_K_Q = (vec_qu_BB*vec_qu_BB');
                mat_N_Q = rho_weight_UL_ite*scal_zu_aux*H_UL_effec*vec_qu_BB';
                % Function to evaluate the optimal matrix of analog combiner/precoder
                [bcd_mat_Qrf_ul] = bcd_analog_matopt(quant_set_ul,mat_Qrf_ul,mat_J_Q,mat_K_Q,mat_N_Q);
                % Update the UL analog combiner
                mat_Qrf_ul = bcd_mat_Qrf_ul;
                
                % DL analog combiner - scalar
                % Search exhaustively in the vector for the one that maximizes the
                % objective function
                obj_func_dl_analog_comb = real((scal_vd_BB'*H_DL_effec'*vec_zd_aux)'*quant_set_dl);
                [~,idx_opt_dl] = max(obj_func_dl_analog_comb);
                % Select the best one and assign it to update the DL analog combiner
                scal_vrf_dl = quant_set_dl(idx_opt_dl);
                
%                 % Update the WMMSE and sum rate for the block
%                 wmse_short_eval;
%                 % Saving the iteration for this block
%                 wmse_eval_iter(iter_BCD,idxBlk_upd) = sum_wmse_total_blk;
%                 rate_eval_iter(iter_BCD,idxBlk_upd) = sum_rate_blk;
%                 % Increase block count
%                 idxBlk_upd = idxBlk_upd + 1;
                
            case 6 % Analog precoders for finite resolution
                %% Analog Precoders Block
                % UL analog precoder - scalar
                % Search exhaustively in the vector for the one that maximizes the
                % objective function
                obj_func_ul_analog_prec = real((scal_zu_aux + delta_unique(iter_Penalty)*scal_lambda_ul(iter_Penalty))'*quant_set_ul);
                [~,idx_opt_ul] = max(obj_func_ul_analog_prec);
                % Select the best one and assign it to update the DL analog combiner
                scal_wrf_ul = quant_set_ul(idx_opt_ul);
                
                % DL analog precoder - matrix
                % Define matrices K_F, and N_F
                mat_K_F = (vec_fbs_BB*vec_fbs_BB');
                mat_N_F = (vec_zd_aux + delta_unique(iter_Penalty)*vec_lambda_dl(:,iter_Penalty))*vec_fbs_BB';
                % Function to evaluate the optimal matrix of analog combiner/precoder
%                 bcd_analog_matopt_asympt_DL;
                [bcd_mat_Frf_dl] = bcd_analog_matopt(quant_set_dl,mat_Frf_dl,eye(par.antBS_Tx),mat_K_F,mat_N_F);
                % Update the UL analog combiner
                mat_Frf_dl = bcd_mat_Frf_dl;
                
%                 % Update the WMMSE and sum rate for the block
%                 wmse_short_eval;
%                 % Saving the iteration for this block
%                 wmse_eval_iter(iter_BCD,idxBlk_upd) = sum_wmse_total_blk;
%                 rate_eval_iter(iter_BCD,idxBlk_upd) = sum_rate_blk;
%                 % Increase block count
%                 idxBlk_upd = idxBlk_upd + 1;
                
            case -1 % Analog combiners for infinite resolution
                %% Analog Combiners Block
                % UL analog combiner - matrix
                % Define matrices J_Q, K_Q, and N_Q
                mat_J_Q = rho_weight_UL_ite*(abs(scal_zu_aux)^2*(H_UL_effec*H_UL_effec') + mat_Psi_ul );
                mat_K_Q = (vec_qu_BB*vec_qu_BB');
                mat_N_Q = rho_weight_UL_ite*scal_zu_aux*H_UL_effec*vec_qu_BB';
                % Function to evaluate the optimal matrix of analog combiner/precoder
                [bcd_mat_Qrf_ul] = bcd_analog_inf_matopt(par.antUE,mat_Qrf_ul,mat_J_Q,mat_K_Q,mat_N_Q);
                % Update the UL analog combiner
                mat_Qrf_ul = bcd_mat_Qrf_ul;
                
                % DL analog combiner - scalar
                % The optimal solution is the complex scalar with modulus equal to
                % par.antBS_Tx
                scal_vrf_dl = (scal_vd_BB'*H_DL_effec'*vec_zd_aux)/...
                    (abs(scal_vd_BB'*H_DL_effec'*vec_zd_aux));%*sqrt(par.antBS_Tx)
                
%                 % Update the WMMSE and sum rate for the block
%                 wmse_short_eval;
%                 % Saving the iteration for this block
%                 wmse_eval_iter(iter_BCD,idxBlk_upd) = sum_wmse_total_blk;
%                 rate_eval_iter(iter_BCD,idxBlk_upd) = sum_rate_blk;
%                 % Increase block count
%                 idxBlk_upd = idxBlk_upd + 1;
                
            case -2 % Analog precoders for infinite resolution
                %% Analog Precoders Block
                % UL analog precoder - scalar
                % The optimal solution is the complex scalar with modulus equal to
                % par.antBS_Tx
                scal_wrf_ul = (scal_zu_aux + delta_unique(iter_Penalty)*scal_lambda_ul(iter_Penalty))/...
                    (abs(scal_zu_aux + delta_unique(iter_Penalty)*scal_lambda_ul(iter_Penalty)));%*sqrt(par.antUE)
                
                % DL analog precoder - matrix
                % Define matrices K_F, and N_F
                mat_K_F = (vec_fbs_BB*vec_fbs_BB');
                mat_N_F = (vec_zd_aux + delta_unique(iter_Penalty)*vec_lambda_dl(:,iter_Penalty))*vec_fbs_BB';
                % Function to evaluate the optimal matrix of analog combiner/precoder
                [bcd_mat_Frf_dl] = bcd_analog_inf_matopt(par.antBS_Tx,mat_Frf_dl,eye(par.antBS_Tx),mat_K_F,mat_N_F);
                % Update the UL analog combiner
                mat_Frf_dl = bcd_mat_Frf_dl;
                
%                 % Update the WMMSE and sum rate for the block
%                 wmse_short_eval;
%                 % Saving the iteration for this block
%                 wmse_eval_iter(iter_BCD,idxBlk_upd) = sum_wmse_total_blk;
%                 rate_eval_iter(iter_BCD,idxBlk_upd) = sum_rate_blk;
%                 % Increase block count
%                 idxBlk_upd = idxBlk_upd + 1;
        end
    end
    %% Evaluate the rates and MSE
    % Updated interference + noise matrix with current beamformer values
    mat_Psi_ul_upd = H_SI*((mat_Frf_dl*vec_fbs_BB)*(mat_Frf_dl*vec_fbs_BB)')*H_SI' +...
        par.noise*eye(par.antBS_Rx);
    
    % Evaluate the UL MSE
    mse_ul_upd = real(abs(1 - (vec_qu_BB'*mat_Qrf_ul')*H_UL_effec*(scal_wrf_ul*scal_wu_BB) )^2 + ...
        (vec_qu_BB'*mat_Qrf_ul')*mat_Psi_ul_upd*(mat_Qrf_ul*vec_qu_BB));
    
    % Evaluate the DL MSE
    mse_dl_upd = real(abs(1 - (scal_vd_BB'*scal_vrf_dl')*H_DL_effec'*(mat_Frf_dl*vec_fbs_BB) )^2 + ...
        (abs(scal_vrf_dl*scal_vd_BB)^2)*scal_Psi_dl);
    
    % UL rate
    sinr_ul = (abs((vec_qu_BB'*mat_Qrf_ul')*H_UL_effec*(scal_wrf_ul*scal_wu_BB) )^2)/...
        real((vec_qu_BB'*mat_Qrf_ul')*mat_Psi_ul_upd*(mat_Qrf_ul*vec_qu_BB));
    rate_ul(iter_BCD) = log2(1 + sinr_ul);
    
    % DL rate
    sinr_dl = ( abs((scal_vd_BB'*scal_vrf_dl')*H_DL_effec'*(mat_Frf_dl*vec_fbs_BB) )^2 )/...
        real( (abs(scal_vrf_dl*scal_vd_BB)^2)*scal_Psi_dl );
    rate_dl(iter_BCD) = log2(1 + sinr_dl);
    
    % Sum rate
    sum_rate(iter_BCD) = rate_ul(iter_BCD) + rate_dl(iter_BCD);
    
    %% Update functions related to BCD convergence
    % UL
    coupling_h_ul_norm_two(iter_BCD+1) = norm(scal_zu_aux - scal_wrf_ul*scal_wu_BB);
    
    % DL
    coupling_h_dl_norm_two(iter_BCD+1) = norm(vec_zd_aux - mat_Frf_dl*vec_fbs_BB);
    
    % Objective function with regularization
    obj_fun_reg_wmmse(iter_BCD) = -log(rho_weight_UL_ite) -log(rho_weight_DL_ite) + ...
        + (rho_weight_UL_ite*mse_ul_upd) + (rho_weight_DL_ite*mse_dl_upd) + ...
        (abs(scal_zu_aux - scal_wrf_ul*scal_wu_BB + delta_unique(iter_Penalty)*scal_lambda_ul(iter_Penalty) )^2)/(2*delta_unique(iter_Penalty)) + ...
        (norm(vec_zd_aux - mat_Frf_dl*vec_fbs_BB + delta_unique(iter_Penalty)*vec_lambda_dl(:,iter_Penalty))^2)/(2*delta_unique(iter_Penalty));
    
    % Objective function without regularization - WMMSE
    obj_fun_wmmse(iter_BCD) = -log(rho_weight_UL_ite) -log(rho_weight_DL_ite) + ...
        + (rho_weight_UL_ite*mse_ul_upd) + (rho_weight_DL_ite*mse_dl_upd);
    
    % Evaluate augmented lagrangian
    aug_Lagrange(iter_BCD) = obj_fun_reg_wmmse(iter_BCD) + real(scal_lambda_ul(iter_Penalty)*...
        (scal_zu_aux - scal_wrf_ul*scal_wu_BB)) + real(vec_lambda_dl(:,iter_Penalty)'*(vec_zd_aux - mat_Frf_dl*vec_fbs_BB) );
    
    %% Evaluate termination criteria
    if iter_BCD > 1
        % Evaluate relative error
        relative_error(iter_BCD-1) = abs(aug_Lagrange(iter_BCD) - aug_Lagrange(iter_BCD-1) )/abs(aug_Lagrange(iter_BCD-1));
        %(abs(sum_rate(iter_BCD) - sum_rate(iter_BCD-1)))/(sum_rate(iter_BCD-1));
        if relative_error(iter_BCD-1) <= epsilon_0
            % Increase the iteration number
            iter_BCD = iter_BCD + 1;
            break;
        end
    end
    % Increase the iteration number
    iter_BCD = iter_BCD + 1;
end
