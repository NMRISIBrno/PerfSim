function config = MakePhantom(config)
% MAKEPHANTOM Load input data, resize them and save to 'tissue' table.
%   Output file: 'tissue.mat' - table containing perfusion and other
%   parameters and binary ROI maps. 
%   Units can be extracted by command: 
%           'tissue.Properties.VariableUnits'
%   Sampling period can be extracted by command:
%           'tissue.Properties.CustomProperties.SamplingPeriod'

    %% Load input data
    wdIn = config.workDir.input;
    [~, ~, ext] = fileparts(config.phantom.roiParam);
    switch ext
        case '.csv'
            if verLessThan('matlab','9.4')
                roiParam =readtable([wdIn, config.phantom.roiParam],'HeaderLines',2);
                headers = readtable([wdIn, config.phantom.roiParam]); % MB - old Matlab 2017b
                roiParam.Properties.VariableUnits = table2cell(headers(1,:));    % define units from 1st line
                roiParam.Properties.VariableNames = headers.Properties.VariableNames;
            else
                roiParam = readtable([wdIn, config.phantom.roiParam],'VariableNamesLine',1,'VariableUnitsLine',2);        
            end
        case '.xlsx'
            roiParam = readtable([wdIn, config.phantom.roiParam],'VariableNamesRange','1:1','VariableUnitsRange','2:2');
        otherwise 
            error('Unknown filetype of input ROI parameter table.')
    end
    roiMap = imread([wdIn, config.phantom.structure]);

    %% Modify input geometry
    numOfRois = roiParam{end,1};
    xSize   = size(roiMap,2);
    ySize   = size(roiMap,1);
    
    if config.phantom.resizeFactor ~= 1
        newMap = imresize(roiMap,config.phantom.resizeFactor,'nearest');
        roiMap = zeros(ySize,xSize);
        yi  = int32(round(ySize-size(newMap,1))/2);
        xi  = int32(round(xSize-size(newMap,2))/2);

        % newmap is placed in the middle with respect to the original map
        roiMap(yi:yi+size(newMap,1)-1, xi:xi+size(newMap,2)-1) = newMap;
    end

    %% Create tissue structure (ROIs are put in cells in table)
    tissue = roiParam;
    tissue.roiMap = cell(numOfRois,1);
    tissue.Properties.VariableUnits{'roiMap'} = '-';
    for roi = 1:numOfRois
        tissue{roi,'roiMap'}{1,1} = roiMap == roi;
    end
    
    % add sampling period to table properties
    if verLessThan('matlab','9.4')
        tissue.Properties.UserData.SamplingPeriod = config.acquisition.timeAxis.Ts;
    else
        tissue = addprop(tissue,'SamplingPeriod','table');
        tissue = addprop(tissue,'SamplingPeriodUnit','table');
        tissue.Properties.CustomProperties.SamplingPeriod = config.acquisition.timeAxis.Ts;
        tissue.Properties.CustomProperties.SamplingPeriodUnit = 's';
    end
    
    %% Save tissue and update config
    config.phantom.xSize = xSize;
    config.phantom.ySize = ySize;
    config.phantom.numOfRois = numOfRois;
    config.phantom.fileName = 'tissue.mat';
    save([config.workDir.output, config.phantom.fileName], 'tissue');   

end
