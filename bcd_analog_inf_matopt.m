function [bcd_mat_analog] = bcd_analog_inf_matopt(num_ant,mat_analog,mat_A,mat_C, mat_B)
%bcd_analog_matopt Function that solves the analog precoder/combiner
%problem using the method proposed by Shi2018, where we update each element
%of the matrix until convergence

% Constant for the BCD
mat_iter = 50;
epsilon_stop = 1e-3;
iter_inner_BCD = 1;
% Define matrix of the analog precoder/combiner
mat_analog_Xtilde = mat_analog;
% Matrix to keep updating while in the BCD loop - Q matrix in Shi2018
mat_loop_Q_upd = mat_A*mat_analog*mat_C;

% Evaluate total objective function
vec_obj_func = zeros(mat_iter,1);
vec_obj_func(1) = real(trace(mat_analog'*mat_A*mat_analog*mat_C) - 2*real(trace(...
    mat_analog'*mat_B)));
relative_diff = ones(mat_iter,1);

% Continue this loop until maximum number
while iter_inner_BCD <= mat_iter
    % Update each element of the matrix
    for idxRow = 1:size(mat_analog,1)
        for idxCol = 1:size(mat_analog,2)
            % Evaluate the scalar that will be used to optimize the matrix
            scal_b = mat_A(idxRow,idxRow)*mat_analog_Xtilde(idxRow,idxCol)*mat_C(idxCol,idxCol) -...
                mat_loop_Q_upd(idxRow,idxCol) + mat_B(idxRow,idxCol);
            % Evaluate the optimal objective function for the scalar
            % problem
            mat_elem = scal_b/(abs(scal_b));%*sqrt(num_ant));
            % Update the element in the Q (Shi2018) matrix
            mat_loop_Q_upd = mat_loop_Q_upd + (mat_elem - mat_analog_Xtilde(idxRow,idxCol))*...
                mat_A(:,idxRow)*mat_C(idxCol,:);
            % Update the desired matrix "mat_analog_bcd"
            mat_analog_Xtilde(idxRow,idxCol) = mat_elem;
        end
    end
    
    % Increase the iteration number
    iter_inner_BCD = iter_inner_BCD + 1;
    % Update the objective function
    vec_obj_func(iter_inner_BCD) = real(trace(mat_analog_Xtilde'*mat_A*mat_analog_Xtilde*mat_C) - 2*real(trace(...
        mat_analog_Xtilde'*mat_B)));
    relative_diff(iter_inner_BCD-1) = abs(vec_obj_func(iter_inner_BCD) - vec_obj_func(iter_inner_BCD-1))/...
        abs(vec_obj_func(iter_inner_BCD-1));
    % Check if the relative difference is lower than predefined (and
    % decreasing) threshold
    if relative_diff(iter_inner_BCD-1) <= epsilon_stop && iter_inner_BCD > 5
        break;
    else
        % Update the threshold
        epsilon_stop = 0.8*epsilon_stop;
    end
end

% Return the variables
bcd_mat_analog = mat_analog_Xtilde;

end