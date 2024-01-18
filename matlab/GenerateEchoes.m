function config = GenerateEchoes(config)
% GENERATEECHOES Generate echo signals based on phantom, coil sensitivities
% and calculated signal intensities.
%   Output file: 'SyntheticEchoes_..._.mat' - echo signals, ksampling noted
%   in filename


    %% Load input data
    wdOut = config.workDir.output;
    wdIn = config.workDir.input;
    load([wdOut, config.outFile.table]);
    load([wdOut, config.acquisition.kSampling.fileName]);
    load([wdIn, config.acquisition.sensitivities]);

    numTEs = length(config.acquisition.TE);
    numOfRois = config.phantom.numOfRois;
    coilElements = size(Sensitivities,3);
    config.acquisition.coilElements = coilElements;
    xSize = config.phantom.xSize;
    ySize = config.phantom.ySize;
    echoes = config.acquisition.timeAxis.N; 
    nSamples = config.acquisition.kSampling.nSamples;
    SD = config.acquisition.SD;
    
    %% Scale sensitivities if normalization constant is given
    if exist('normConst','var')
        Sensitivities = Sensitivities ./ normConst;
    end
    
    %% Resize sensitivities to the same size (FOV) as phantom
    if config.phantom.resizeFactor ~= 1
        newMap = imresize(Sensitivities,config.phantom.resizeFactor,'nearest');
        Sensitivities = zeros(ySize,xSize,coilElements);
        yi  = int32(round(ySize-size(newMap,1))/2);
        xi  = int32(round(xSize-size(newMap,2))/2);

        % newmap is placed in the middle with respect to the original map
        Sensitivities(yi:yi+size(newMap,1)-1, xi:xi+size(newMap,2)-1, :) = newMap;
    end
    
    %% Generate echo signals
    echoSignals = zeros(nSamples,echoes,numTEs,coilElements);
    
    k = complex(ky,kx); % samplesPerEcho x numberOfEchoesPerTE  ...ugly trick for NUFFT and dtft
    om = [real(k(:))*ySize, imag(k(:))*xSize, zeros(size(k(:)))]'; % totalNumberOfSampledKSpaceSamples x 3 (kx,ky,kz directions)
    
    fprintf('Generating data:\n')
    
    switch config.acquisition.kSampling.method
        case 'radial' 
            for TEind = 1:numTEs
                obj = nufft_3d(om,[ySize,xSize,1]');
                for coil = 1:coilElements
                    for roi = 1:numOfRois
                        fprintf('TEind: %d/%d, Coil: %d/%d, ROI: %d/%d\n', TEind,numTEs, coil,coilElements, roi,numOfRois)
                        dataOut = obj.fNUFT(tissue{roi,'roiMap'}{1,1}.*Sensitivities(:,:,coil)); % totalNumberOfSampledKSpaceSamples x 1
                        dataOut = reshape(dataOut,size(k)); % samplesPerEcho x numberOfEchoesPerTE
                        SI = reshape(tissue{roi,'SI'}{1}(TEind,:), [1 echoes]); % 1 x numberOfEchoesPerTE
                        SI = repmat(SI,size(dataOut,1),1); % samplesPerEcho x numberOfEchoesPerTE
                        dataOut = dataOut.*SI;
                        echoSignals(:,:,TEind,coil) = echoSignals(:,:,TEind,coil) + dataOut;
                    end
                end
            end
        case 'cartesian'
            phEncSteps = config.acquisition.kSampling.phEncSteps;
            for TEind = 1:numTEs
                for coil = 1:coilElements
                    for roi = 1:numOfRois
                        fprintf('TEind: %d/%d, Coil: %d/%d, ROI: %d/%d\n', TEind,numTEs, coil,coilElements, roi,numOfRois)
                        image = imresize(tissue{roi,'roiMap'}{1,1}.*Sensitivities(:,:,coil),[phEncSteps nSamples],'Method','bilinear');
                        dataOut = fftshift(fft2(fftshift(image)/numel(image)))'; 
                        dataOut = repmat(dataOut,1,ceil(echoes/phEncSteps)); % copy images next to each other
                        dataOut = dataOut(:,1:echoes);
                        SI = reshape(tissue{roi,'SI'}{1}(TEind,:), [1 echoes]);
                        SI = repmat(SI,size(dataOut,1),1);
                        dataOut = dataOut.*SI; % echoes are columns
                        echoSignals(:,:,TEind,coil) = echoSignals(:,:,TEind,coil) + dataOut;
                    end
                end
            end
    end
    
    % add zero-mean Gaussian noise
    echoSignals = echoSignals + SD*randn(size(echoSignals)) + 1i*SD*randn(size(echoSignals));
    
    %% Save echo data
    if ~isfield(config.outFile,'echoes')
        config.outFile.echoes = ['SyntheticEchoes_', config.acquisition.kSampling.method, '.mat'];
    end
    save([wdOut, config.outFile.echoes], 'echoSignals')
    
    %% Save output config as .json
    txt = jsonencode(config);
    fileName = 'config.json';
    fid = fopen([wdOut, fileName],'w');
    fprintf(fid, "%s", txt);
    fclose(fid);
    
end