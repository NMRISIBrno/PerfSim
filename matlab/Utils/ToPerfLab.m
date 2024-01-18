function [data1, data2, info] = ToPerfLab(path, precontrast)

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

    % Prealocate data1
    data1 = cell(1,length(pars));
    info.acq.kM0Range = size(data1); 

    % Process precontrast data
    for prec = 1:length(pars)
        % Load echoes
        load([path,'/Data1/SyntheticEchoes_cartesian_',config.acquisition.calibration.model,'_',num2str(pars(prec))])
        
        % Extract some parameters from config
        phEncSteps = config.acquisition.kSampling.phEncSteps;
        frames = config.acquisition.kSampling.repetitions;
        nSamples = config.acquisition.kSampling.nSamples;
        coilElements = config.acquisition.coilElements;
        numOfEchoes = length(config.acquisition.TE);

        % Prepare data and vars
        totalSamples = phEncSteps*frames;
        echoSignals = echoSignals(:,1:totalSamples,:,:);
        echoSignals = reshape(echoSignals, nSamples, phEncSteps, frames, numOfEchoes, coilElements);
        
        imdata = zeros(phEncSteps, nSamples, frames, numOfEchoes, coilElements);
        data = cell(1,frames*numOfEchoes);

        % Plot reconstructed data
        figure(1)
        for i = 1:frames
            for j = 1:numOfEchoes
                for k = 1:coilElements
                    imdata(:,:,i,j,k) = abs(fftshift(ifft2(fftshift(echoSignals(:,:,i,j,k)')))).^2; % squared abs image
                end
                data{1,(i-1)*numOfEchoes+j} = sqrt(sum(imdata(:,:,i,j,:),5));
                imagesc(data{1,(i-1)*numOfEchoes+j});
                title(['Image: ', num2str((i-1)*numOfEchoes+j)])
                pause(0.01);
            end
        end

        % Save to data1
        data1{1, prec} = data;
    end

    % Save parameters for info structure
    if strcmp(config.acquisition.calibration.model,'vFA')
        TR = repmat(config.acquisition.TR*1000,1,length(pars));      % from [s] to [ms]
        FA = pars;   
    elseif strcmp(config.acquisition.calibration.model,'vTR')
        TR = pars;      
        FA = repmat(config.acquisition.FA./pi*180,1,length(pars));   % from [rad] to [deg]
    else
        TR = repmat(config.acquisition.TR*1000,1,length(pars));      % from [s] to [ms]
        FA = repmat(config.acquisition.FA./pi*180,1,length(pars));   % from [rad] to [deg]
    end
    TE = repmat(config.acquisition.TE*1000,1,length(pars));      % from [s] to [ms]
else
    % Otherwise create empty data1
    data1 = {};
end


%% Get dynamic data and create data2 structure
configPath = [path,'/Data2/config.json'];

% Load config
fid = fopen(configPath);    % opening the file
raw = fread(fid,inf);       % reading the contents
str = char(raw');           % transformation
fclose(fid);

config = jsondecode(str);

% Load echoes
load([path,'/Data2/SyntheticEchoes_cartesian'])

% Extract some parameters from config
phEncSteps = config.acquisition.kSampling.phEncSteps;
frames = config.acquisition.kSampling.repetitions;
nSamples = config.acquisition.kSampling.nSamples;
coilElements = config.acquisition.coilElements;
numOfEchoes = length(config.acquisition.TE);

% Prepare data and vars
totalSamples = phEncSteps*frames;
echoSignals = echoSignals(:,1:totalSamples,:,:);
echoSignals = reshape(echoSignals, nSamples, phEncSteps, frames, numOfEchoes, coilElements);

imdata = zeros(phEncSteps, nSamples, frames, numOfEchoes, coilElements);
data2 = cell(1,frames*numOfEchoes);

% Plot reconstructed data
figure(2)
for i = 1:frames
    for j = 1:numOfEchoes
        for k = 1:coilElements
            imdata(:,:,i,j,k) = abs(fftshift(ifft2(fftshift(echoSignals(:,:,i,j,k)')))).^2; % squared abs image
        end
        data2{1,(i-1)*numOfEchoes+j} = sqrt(sum(imdata(:,:,i,j,:),5));
        imagesc(data2{1,(i-1)*numOfEchoes+j});
        title(['Image: ', num2str((i-1)*numOfEchoes+j)])
        pause(0.01);
    end
end

% Save parameters for info structure
TR = [TR, config.acquisition.TR*1000];      % from [s] to [ms]
TE = [TE, config.acquisition.TE*1000];      % from [s] to [ms]
FA = [FA, config.acquisition.FA./pi*180];   % from [rad] to [deg]
        

%% Create remaining info structure
info.acq.TR = TR;
info.acq.TE = TE;
info.acq.FA = FA;
info.acq.Ts = phEncSteps * config.acquisition.timeAxis.Ts;

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