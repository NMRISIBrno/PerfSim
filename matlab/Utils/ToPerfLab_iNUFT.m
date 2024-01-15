function [data1, data2, info] = ToPerfLab_iNUFT(path, precontrast, nRadial)

% Get folder names
folders = dir(path);

% Prealocate info parameters
TR = [];
TE = [];
FA = [];

%% Get precontrast data and create data1 structure
if precontrast
    % Load config
    configPath = [path,'/Data1/config.json'];

    fid = fopen(configPath);    % opening the file
    raw = fread(fid,inf);       % reading the contents
    str = char(raw');           % transformation
    fclose(fid);

    config = jsondecode(str);
        
    pars = config.acquisition.calibration.pars';
    
    % Load k-space trajectory
    load([path,'/Data1/',config.acquisition.kSampling.fileName])

    % Prealocate data1
    data1 = cell(1,length(pars));
    info.acq.kM0Range = size(data1); 

    % Process precontrast data
    for prec = 1:length(pars)
        % Load echoes
        load([path,'/Data1/SyntheticEchoes_radial_',config.acquisition.calibration.model,'_',num2str(pars(prec))])
 
        % Extract some parameters from config
        projections = config.acquisition.kSampling.projections;
        frames = floor(projections/nRadial);

        for f=1:frames
           kxf = kx(:,((f-1)*nRadial+1):f*nRadial);
           kyf = ky(:,((f-1)*nRadial+1):f*nRadial);

           om = [kyf(:)*1024, kxf(:)*1024, zeros(size(kxf(:)))]';
           obj = nufft_3d(om, [128,128,1]');

           data1{prec}{f} = abs(obj.iNUFT(echoSignals(:,((f-1)*nRadial+1):f*nRadial),10,0,0));
           
           figure(1);
           imagesc(data1{prec}{f})
           pause(0.01);
        end

        % Save parameters for info structure
        TR = [TR, config.acquisition.TR*1000];      % from [s] to [ms]
        TE = [TE, config.acquisition.TE*1000];      % from [s] to [ms]
        FA = [FA, config.acquisition.FA./pi*180];   % from [rad] to [deg]
    end
    
else
    % Otherwise create empty data1
    data1 = {};
end


%% Get dynamic data and create data2 structure
configPath = [path,'Data2/config.json'];

% Load config
fid = fopen(configPath);    % opening the file
raw = fread(fid,inf);       % reading the contents
str = char(raw');           % transformation
fclose(fid);

config = jsondecode(str);

% Load echoes and k-sampling
load([path,'/Data2/SyntheticEchoes_radial'])
load([path,'/Data2/',config.acquisition.kSampling.fileName])

% Extract some parameters from config
projections = config.acquisition.kSampling.projections;
frames = floor(projections/nRadial);

for f=1:frames
   kxf = kx(:,((f-1)*nRadial+1):f*nRadial);
   kyf = ky(:,((f-1)*nRadial+1):f*nRadial);

   om = [kyf(:)*1024, kxf(:)*1024, zeros(size(kxf(:)))]';
   obj = nufft_3d(om, [128,128,1]');

   data2{f} = abs(obj.iNUFT(echoSignals(:,((f-1)*nRadial+1):f*nRadial),10,0,0));

   figure(2);
   imagesc(data2{f})
   pause(0.01);
end

% Save parameters for info structure
TR = [TR, config.acquisition.TR*1000];      % from [s] to [ms]
TE = [TE, config.acquisition.TE*1000];      % from [s] to [ms]
FA = [FA, config.acquisition.FA./pi*180];   % from [rad] to [deg]
        

%% Create remaining info structure
info.acq.TR = TR;
info.acq.TE = TE;
info.acq.FA = FA;

info.section = 1;
info.sections = size(data2,1);
info.codewords = {'inp'};
info.source = {configPath};
info.modality = 'mri';

info.author = 'Simulator2.0';
info.experiment = 'PhantomRat';
info.notes = {''}; % name of PK model not available
% info.acq.kM0Range = [1 length(info.acq.FA)-1];
[filename, info] = name_file(info);

%% Save file
save([path, filename],'data1','data2','info')

end