% Needs to know the number of files and the number of channels
% Perhaps load the first file to get the number of channels and 
% channel names?

alpha = 0.05; % the significance level of the connectivity between sites




% This part initializes the output and writes the header. 
% Needs an EEG file loaded with the number of channels and channel names
OutputMatrix = cell(NumFile+1, NumChan.^2+1);
OutputMatrix{1,1} = 'Filename';
count=1;
for i=1:NumChan
    for j=1:NumChan
        count=count+1;
        OutputMatrix{1,count} = strcat(EEG.chanlocs(i).labels, '_to_', EEG.chanlocs(j).labels);
    end
end

% This part generates the actual results for each EEG file
for i=1:NumFile
    % Load the EEG structure for file i
    OutputMatrix{i+1, 1} = EEG.filename;
    Statistic = EEG.PSGC.Statistic;
    theM = EEG.PSGC.theM;
    
    M=repmat(theM,NumChan,1);
    critical_value = zeros(NumChan);
    for j=1:NumChan
        for k=1:NumChan
            critical_value(j,k) = chi2inv(1-alpha, M(j,k));
        end
    end
    Connected = Statistic > critical_value;
    count=1;
    for j=1:NumChan
        for k=1:NumChan
	    count=count+1;
            OutputMatrix{i+1, count} = Connected(j,k);
        end
    end
end
            
    
