clear variables
delete(gcp('nocreate'))

c = parcluster('tetralith');
c.AdditionalProperties.AccountName = 'snic2021-1-2';
c.AdditionalProperties.WallTime = '3-00:00:00';
c.AdditionalProperties.EmailAddress = 'x_marma@tetra.nsc.liu.se';
c.saveProfile;

addpath(genpath('/home/x_marma/Ua/UaSource/'))
addpath(genpath('/home/x_marma/Ua/PICO_Ua/'))
addpath(genpath('/home/x_marma/Ua/scripts/'))

numOutputs = 1;
numWorkers = 31; % n = numWorkers+1

cd /home/x_marma/Ua/runs/jutul_slope_proj_sl12p7_calv/
c.AdditionalProperties.AdditionalSubmitArgs = '--job-name=sl12p7';
j = batch(c, @Ua, numOutputs, 'Pool', numWorkers, 'CurrentDirectory', '.');

% j.diary will show the screen outputs

