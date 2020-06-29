%% Parameters
% JointBF                         = 1;      % LOWRES_PHASE or INFINITE_PHASE
% AntennaAlloc                    = 1;      % Antenna Splitting
% lamdamu                         = [5];    % Total number of UL users
% lambdadl                        = 6;      % Total number of DL users
% antBS                           = 64;     % Number of antennas at the BS
% antBS_RF                        = 2       % Number of RF chains at the BS
% antUE                           = 4;      % Number of antennas at the BS
% antUE_RF                        = 1;      % Number of RF chains at the UE
% beta                            = 110;    % SI cancellation [-dB]
% seedMC                          = 1:100;  % Monte Carlo iteration

% Size of matrices
% H_UL = [ h_1[Mbs X Mue] ... h_lu[Mbs X Mue] ] => H_UL [Mbs X (Mue*lu)]
% H_DL = [ h_1[Mbs X Mue] ... h_ld[Mbs X Mue] ] => H_DL [Mbs X (Mue*ld)]
% H_mm = [ h_{u1d1}[Mue X Mue] ... h_{u1 ld}[Mue X Mue] ] => H_mm [(Mue*lu) X (Mue*ld)]