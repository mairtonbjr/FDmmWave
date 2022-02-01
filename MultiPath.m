function [ H_effect ] = MultiPath(L_multi,g,ant_Tx,ant_Rx)
%MultiPath Creates the multipath channel from the array steerging
%vectors

% Number of users
N_user = size(g,1);

% Generate the TX and RX angles for each path
angle_tx = 2*pi*rand(N_user,L_multi) - pi;
angle_rx = 2*pi*rand(N_user,L_multi) - pi;

% Define zero vector with this size
H_effect = zeros(ant_Tx,1);

% For to generate the channel for each user
for idxUser = 1:N_user
    % If the channel is already a matrix, then H_effect should be a cell
    % For now, H_effect is a vector, and since we have only 1 UL/DL user,
    % H_effect is a vector
    for idxLpath = 1:L_multi
        % Vector to span exponential array
        vec_span_complex_tx = exp(1j*(0:ant_Tx-1)*(pi*sin(angle_tx(idxLpath)) ) );
        vec_span_complex_rx = exp(1j*(0:ant_Rx-1)*(pi*sin(angle_rx(idxLpath)) ) );
        % Generate the TX and RX array responses
        arra_resp_tx = sqrt(1/ant_Tx)*vec_span_complex_tx;
        arra_resp_rx = sqrt(1/ant_Rx)*vec_span_complex_rx;
        
        % Generate the channel
        H_effect = H_effect + g(idxLpath)*arra_resp_rx*arra_resp_tx';
    end
end

% Final product that is independent of the sums above
H_effect = sqrt(ant_Tx*ant_Rx/L_multi)*H_effect;
end