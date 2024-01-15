function config = PKmodeling(config)
% PKMODELING Load input data and calculate SI based on PK models given for
% particular ROIs. Add all curves to the 'tissue' table.
%   Output file: 'tissue.mat' - table with input parameters and calculated
%   signals


    %% Load input data
    wdOut = config.workDir.output;
    load([wdOut, config.phantom.fileName]);
    load([wdOut, config.phantom.AIF.fileName])

    numOfRois = config.phantom.numOfRois;
    N = config.acquisition.timeAxis.N;
    numTEs = length(config.acquisition.TE);

    %% PK modeling
    % create time axis
    t  = (0:(N-1))*config.acquisition.timeAxis.Ts/60; % [min]
    allCurves = [];
    
    for roi = 1:numOfRois
       
        % get PK model for specific ROI
        fh = str2func(tissue{roi,'PK_model'}{1,1});
        curves = fh(config, aif, t, tissue(roi,:)); 
        
        % get SI for all TEs
        SIs = zeros(numTEs, N);
        fh = str2func(config.acquisition.model);
        for TEind = 1:numTEs
            [SI, SITE0] = fh(config, curves.R1{1,1}, curves.R2{1,1}, TEind);
            SIs(TEind,:) = SI;
        end
        
        curves.SITE0 = {SITE0};
        curves.SI = {SIs};

        allCurves = [allCurves; curves];
    end
    
    % Add curves to tissue table
    tissue = [tissue allCurves];
    tissue.Properties.VariableUnits{'SITE0'} = '-';    
    tissue.Properties.VariableUnits{'SI'} = '-';
    
    %% Update config and save table with all parameters and curves
    if ~isfield(config.outFile,'table')
        config.outFile.table = 'tissue.mat';
    end
    save([wdOut, config.outFile.table], 'tissue');

end