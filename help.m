%% Parameters
% JointBF                         = LOWRES_PHASE;      % LOWRES_PHASE or INFINITE_PHASE
% AntennaAlloc                    = SplittingFD;       % Antenna Splitting between UL and DL
% lamdamu                         = 1;                 % Total number of UL users
% lambdadl                        = 1;                 % Total number of DL users
% antBS                           = 64;                % Number of antennas at the BS
% antBS_RF                        = 2                  % Number of RF chains at the BS
% antUE                           = 1;                 % Number of antennas at the UE
% antUE_RF                        = 1;                 % Number of RF chains at the UE
% beta                            = 110;               % SI cancellation [-dB]
% seedMC                          = 1:100;             % Monte Carlo iteration


%% Define important parameters within the structure par
% Assign the antennas and RF chains
par.antBS = 64;
par.antBS_Tx = par.antBS/2;
par.antBS_Rx = par.antBS/2;
par.antBS_RF = 2;
par.antUE = 1;
par.antUE_RF = 1;
% SI cancellation in linear scale
par.beta = sqrt(db2lin(-110));
% Seed
par.seedMC = 1;
% Number of UL and DL users
par.lambdaul = 1;
par.lambdadl = 1;
% Set the general parameters
par.JointBF = 'LOWRES_PHASE'; % LOWRES_PHASE or INFINITE_PHASE
% Antenna Splitting between UL and DL. If 64 antennas in total, use 32 for
% UL and 32 for DL
par.AntAlloc = 'SplittingFD';
% Phase shifter resolution
par.bit_res = 6; % 1, 3, 6, or inf

% Transmission specific parameters 
par.pmaxUL = db2lin(23-30); % maximum uplink power
par.pmaxDL = db2lin(24-30); % maximum downlink power
% Asymptotic analysis is not present in this paper
par.asympt = 0;
% No frequency selective fading
par.FreqFad = 0;
% Channel bandwidth
par.chbw = 720e3;%180e3; %Bandwidth of each PRB - KHz 60*12 for mmWave
% Noise power
par.noise = (10^(-17.4))*par.chbw*db2lin(13); % High noise figure


%% Generate the channels

% Size of matrices
% H_UL = [ h_1[antBS X antUE] ... h_lu[antBS X antUE] ] => H_UL [antBS X (antUE*lu)]
% H_DL = [ h_1[antBS X antUE] ... h_ld[antBS X antUE] ] => H_DL [antBS X (antUE*ld)]
% H_mm = [ h_{u1d1}[antUE X antUE] ... h_{u1 ld}[antUE X antUE] ] => H_mm [(antUE*lu) X (antUE*ld)]
% H_SI [antBS_Rx X antBS_Tx]

% Type of pathloss channel used in the paper is based on the paper:
% M. R. Akdeniz et al., "Millimeter Wave Channel Modeling and
% Cellular Capacity Evaluation," in IEEE Journal on Selected Areas
% in Communications, vol. 32, no. 6, pp. 1164-1179, June 2014.
par.channel = 'NewYork14'; % '3GPP', 'WINNER_O2Ia', 'WINNER_UMi' (usual), '3GPP_Pico', 'Uniform', 'NewYork14'

% Number of paths for the multi-path array steering
par.L_multi = 3; % 3, 6, 12
% Rician power factor in dB
par.K_Rice = db2lin(50); % 0, 30, 50
par.SI_model = 'Ricean';

% Necessary input to generate channels
% This is the pathloss in linear scale. You should choose the pathloss model 
% to the scenario of your application. I am not including the
% implementation here, so the pathloss here is simply 1.
gUL = ones(par.lambdaul,1); 
gDL = ones(par.lambdadl,1); 
% Include the L multipaths - complex gains with unitary variance
gUL = sqrt(repmat(gUL,1,par.L_multi)/2).*complex(randn(par.lambdaul,par.L_multi),...
    randn(par.lambdaul,par.L_multi));
gDL = sqrt(repmat(gDL,1,par.L_multi)/2).*complex(randn(par.lambdadl,par.L_multi),...
    randn(par.lambdadl,par.L_multi));

% UL and DL using the multipath model
H_UL_effec = MultiPath(par.L_multi,gUL,par.antBS_Tx,par.antUE);
H_DL_effec = MultiPath(par.L_multi,gDL,par.antUE,par.antBS_Tx);

% SI
cov_mat = chol( 1/(1+par.K_Rice)*(eye(par.antBS_Rx,par.antBS_Tx)) );
H_SI = (1/sqrt(2))*(complex(randn(par.antBS_Rx,par.antBS_Tx),randn(par.antBS_Rx,par.antBS_Tx))*cov_mat + ...
    sqrt((par.K_Rice)/(1+par.K_Rice))*ones(par.antBS_Rx,par.antBS_Tx));
H_SI = H_SI*par.beta;%/sqrt(par.antBS_Rx*par.antBS_Tx)

%% Generate the folder to store the results

% Define a name for the MC iteration. 
[par] = setFileName(par);

%% Run the simulations

[obj_fun_reg_wmmse,vec_coupling,sum_rate_PDD,sum_rate,SpEff_final,beam_var] = wmmse_hybrid_precod(par,H_UL_effec,H_DL_effec,H_SI);