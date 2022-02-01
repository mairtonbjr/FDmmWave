function [par] = setFileName(par)

% Results directory
folder = './results_fd/';

if par.FreqFad
    fading = strcat('FreqFad_',par.channel,'_');
else
    fading = strcat(par.channel,'_');
end

if ~isinf(par.bit_res)
    bit_res_str = strcat('_BitRes',num2str(par.bit_res),'_');
else
    bit_res_str = '_';
end

if strcmp(par.SI_model,'Ricean')
    if par.K_Rice ~= 1
        SI_channel = strcat(num2str(lin2db(par.beta^2)),'dB_K',num2str(lin2db(par.K_Rice)),'dB');
    else
        SI_channel = strcat(num2str(lin2db(par.beta^2)),'dB');
    end
end

% Folder
par.save_path = strcat(folder,fading,SI_channel,bit_res_str,...
    par.AntAlloc,'_',par.JointBF,'_PwrDL',num2str(lin2dbm(par.pmaxDL)),'dBm_',num2str(par.lambdaul),'UL_',...
    num2str(par.lambdadl),'DL_',strcat(num2str(par.antBS),'-',num2str(par.antBS_RF),'MBS_'),...
    strcat(num2str(par.antUE),'-',strcat(num2str(par.antUE_RF),'MUE/')));

% File name
par.par_save_name = strcat(par.save_path,'par','.mat');
par.sta_save_name = strcat(par.save_path,'sta_glob');
if ~exist(strcat(par.sta_save_name,'_',num2str(par.seedMC),'.mat'), 'file')
    if ~exist(par.save_path, 'dir')
        mkdir(par.save_path);
    end
else
    disp('This simulation has already been done.');
end

end
