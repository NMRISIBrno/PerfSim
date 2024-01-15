function outputData = GetParameterMaps(tissue, outputVars, outputSize, outputTs)
% GETPARAMETERMAPS Create perfusion parameter maps from input data.
%   Input:  tissue     - table with rois and parameters
%           outputVars - cell with names of output parameters
%               (calculated parameters: vp, Ktrans, kep, PS)
%               'scalar'  - Fp, E, ve, Tc, vp, Ktrans, kep, PS
%               'dynamic' - TRF, ct, R1, R2, SITE0, SI
%               'all'         - scalar + dynamic
%           outputSize - size of maps, e.g. [128 128]
%           outputTs   - new sampling period for vector parameters [s]
%   Output: outputData - table with generated perfusion parameter maps


if strcmp(outputVars, 'scalar')
    outputVars = {'Fp', 'E', 've', 'Tc', 'vp', 'Ktrans', 'kep', 'PS'};
elseif strcmp(outputVars, 'dynamic')
    outputVars = {'TRF', 'ct', 'R1', 'R2', 'SITE0', 'SI'};
elseif strcmp(outputVars, 'all')
    outputVars = {'Fp', 'E', 've', 'Tc', 'vp', 'Ktrans', 'kep', 'PS', ...
        'TRF', 'ct', 'R1', 'R2', 'SITE0', 'SI'};
end

numOfRois = size(tissue,1);
numOfVars = length(outputVars);
xSize = size(tissue{1,'roiMap'}{1,1},2);
ySize = size(tissue{1,'roiMap'}{1,1},1);
varTypes = repmat({'cell'},1,numOfVars);
Ts = tissue.Properties.CustomProperties.SamplingPeriod;

% set output size (if not given) as the input
if ~exist('outputSize','var')
    outputSize = [ySize, xSize];
end

% set output Ts (if not given) as the input
if ~exist('outputTs','var')
    outputTs = Ts;
end

% resampling parameters
N = length(tissue{1,'TRF'}{1,1});
outputN = floor(Ts*N/outputTs);
x = 0:Ts:(N-1)*Ts;                      % time axis, sample points
xq = 0:outputTs:(outputN-1)*outputTs;   % time axis, query points

% create output table
outputData = table('Size', [1,numOfVars], 'VariableTypes', varTypes, 'VariableNames', outputVars);


%% Resize ROIs
rois = zeros([outputSize,numOfRois]);
phantom = zeros(ySize,xSize);

% create phantom
for roi = 1:numOfRois
    mask = tissue{roi,'roiMap'}{1,1};
    phantom(mask) = roi;
end

% resize it
phantom = imresize(phantom,outputSize,'nearest');

% split to logical 3D matrix
for roi = 1:numOfRois
    mask = phantom == roi;
    rois(:,:,roi) = mask;
end
rois = logical(rois);


%% Calculate parameter maps
for var = 1:numOfVars 
    
    % Check if variable exists in table
    if any(ismember(tissue.Properties.VariableNames, outputVars{var}))
    
        % For scalar parameters
        if isnumeric(tissue{1,outputVars{var}})
            
            map = zeros(outputSize);
            % go through all rois and set parameter in map
            for roi = 1:numOfRois
                map(rois(:,:,roi)) = tissue{roi,outputVars{var}};
            end
            % save map with units 
            outputData{1,outputVars{var}}{1,1} = map;
            outputData.Properties.VariableUnits{outputVars{var}} = tissue.Properties.VariableUnits{outputVars{var}};


        % For vector parameters
        elseif iscell(tissue{1,outputVars{var}}) && isvector(tissue{1,outputVars{var}}{1,1}) && isnumeric(tissue{1,outputVars{var}}{1,1})

            maps = zeros([outputSize,outputN]);

            % go through all rois and resample parameter
            for roi = 1:numOfRois
                
                % extract vector parameter for current ROI
                vec = tissue{roi,outputVars{var}}{1,1};
                
                % resample vector
                vec = interp1(x,vec,xq,'linear');
                
                % set values in dynamic sequence
                mask3d = repmat(rois(:,:,roi),1,1,outputN);
                vec = repmat(vec,sum(rois(:,:,roi),'all'),1);
                maps(mask3d) = vec;
            end  

            outputData{1,outputVars{var}}{1,1} = maps;
            outputData.Properties.VariableUnits{outputVars{var}} = tissue.Properties.VariableUnits{outputVars{var}};

        % For resampled roi map    
        elseif strcmp(outputVars{var},'roiMap')
            outputData{1,outputVars{var}}{1,1} = phantom;
            
            
        % Parameter is in table, but it's not a scalar/vector
        else
            error([outputVars{var}, ': Cannot create map from this datatype. Choose scalar/vector parameter.'])
        end

        
    % If variable isn't in table, calcute it
    else
        switch outputVars{var}
            case 'vp'
                vp = zeros(ySize,xSize);
                % go through all rois and set parameter in map
                for roi = 1:numOfRois
                    mask = tissue{roi,'roiMap'}{1,1};
                    vp(mask) = tissue{roi,'Tc'} .* tissue{roi,'Fp'};
                end
                % resize map
                vp = imresize(vp,outputSize,'nearest');
                % save map with units
                outputData{1,'vp'}{1,1} = vp;
                if strcmp(tissue.Properties.VariableUnits{'Tc'},'min') && strcmp(tissue.Properties.VariableUnits{'Fp'},'ml/min/ml')
                    outputData.Properties.VariableUnits{'vp'} = 'ml/ml';
                else
                    outputData.Properties.VariableUnits{'vp'} = 'unknown';
                end

            case 'Ktrans'
                Ktrans = zeros(ySize,xSize);
                % go through all rois and set parameter in map
                for roi = 1:numOfRois
                    mask = tissue{roi,'roiMap'}{1,1};
                    Ktrans(mask) = tissue{roi,'Fp'} .* tissue{roi,'E'};
                end
                % resize map
                Ktrans = imresize(Ktrans,outputSize,'nearest');
                % save map with units
                outputData{1,'Ktrans'}{1,1} = Ktrans;
                if strcmp(tissue.Properties.VariableUnits{'Fp'},'ml/min/ml') && strcmp(tissue.Properties.VariableUnits{'E'},'-')
                    outputData.Properties.VariableUnits{'Ktrans'} = '1/min';
                else
                    outputData.Properties.VariableUnits{'Ktrans'} = 'unknown';
                end

            case 'kep'
                kep = zeros(ySize,xSize);
                % go through all rois and set parameter in map
                for roi = 1:numOfRois
                    mask = tissue{roi,'roiMap'}{1,1};
                    kep(mask) = tissue{roi,'E'} .* tissue{roi,'Fp'} ./ (tissue{roi,'ve'} + 1e-20); % to avoid NaNs
                end
                % resize map
                kep = imresize(kep,outputSize,'nearest');
                % save map with units
                outputData{1,'kep'}{1,1} = kep;
                if strcmp(tissue.Properties.VariableUnits{'E'},'-') && strcmp(tissue.Properties.VariableUnits{'Fp'},'ml/min/ml') && strcmp(tissue.Properties.VariableUnits{'ve'},'ml/ml')
                    outputData.Properties.VariableUnits{'kep'} = '1/min';
                else
                    outputData.Properties.VariableUnits{'kep'} = 'unknown';
                end

            case 'PS'
                PS = zeros(ySize,xSize);
                % go through all rois and set parameter in map
                for roi = 1:numOfRois
                    mask = tissue{roi,'roiMap'}{1,1};
                    % equation based on model - plug-flow/compartment
                    switch tissue{roi,'PK_model'}{1}
                        case {'ATHv01','aprox_aaTH_DCATH_v5_DceDsc','TH_trunc_FT'}
                            PS(mask) = -log(1-tissue{roi,'E'}) .* tissue{roi,'Fp'};
                        case {'TwoCX','TwoCU_FT'}
                            PS(mask) = tissue{roi,'E'} .* tissue{roi,'Fp'} ./ (1-tissue{roi,'E'});
                        otherwise
                            error('Undefined PS equation for used PK model.')
                    end
                end
                % resize map
                PS = imresize(PS,outputSize,'nearest');
                % save map with units
                outputData{1,'PS'}{1,1} = PS;
                if strcmp(tissue.Properties.VariableUnits{'E'},'-') && strcmp(tissue.Properties.VariableUnits{'Fp'},'ml/min/ml')
                    outputData.Properties.VariableUnits{'PS'} = 'ml/min/ml';
                else
                    outputData.Properties.VariableUnits{'PS'} = 'unknown';
                end
            
            % No option for the parameter
            otherwise
                error([outputVars{var}, ': Unknown parameter.'])
                
        end
    end
end


end