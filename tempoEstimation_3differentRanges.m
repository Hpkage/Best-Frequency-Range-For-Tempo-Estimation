%% Load Data

X = loadMultipleFiles();  % Feature vector

nComp = 1; % How many components to return
nReg = 7; % Regularization parameter - probably don't need to change ever 

%% Spatially Filter & Extract Useful Data

% Finds the most useful signal information in the dataset
% We only care about dataOut (data of the EEG that was filtered)
% Out of the 125 electric cusps used to capture brain waves, 
[dataOut,W,A,Rxx,Ryy,Rxy,dGen] = rcaRun125_parpoolAlready2021(...
    X, nReg, nComp, [], [], 0, []);

%% Convert to 2D Matrix & Plot

for i = 1:10
    % squeeze converts the 3D matrix into a 2D matrix
    % fft converts discrete signals from the time domain into a frequency domain
    dataFFT{i} = abs(fft(squeeze(dataOut{i})));
    dataMean{i} = abs(fft(mean(squeeze(dataOut{i}),2)));
    fAx{i} = computeFFTFrequencyAxis(size(dataFFT{i},1), 125);
    
    % Plot every participants' EEG data
    % figure();
    % plot(fAx{i}, dataFFT{i});
    % hold on    % Combines the plots
    % plot(fAx{i}, dataMean{i}, 'k', 'linewidth', 2);   % Add mean into the plot (black peaks)
    % xlim([0,15]);   % Limit the frequency from 0 to 15 Hz 
    % ylim([0,10*10^4]);
    % set(gca,'YTickLabel',[]);
    % title('EEG Data: Song', i);
    % xlabel('Frequency (Hz)');
    % ylabel('Magnitude (A.U.)');
end

%% Get Largest Peak Between 3 to (6+k) Hertz

% Using the largest peak values for ranges between 3-7 Hz, 3-8 Hz, & 3-9 Hz
% to analyze which range produces the best result

for k = 1
    for i = 9
        % finds index of values between 3 to (6+k) Hz
        idxUse = find(fAx{i} >= 0 & fAx{i} <= 15);
    
        % subset data values between 3 to (6+k) Hz
        fftUse = dataFFT{i}(idxUse,:);
    
        % create matrix of largest peak observed between 3 to (6+k) Hz
        for j = 1:20
            peakIdx = find(fftUse(:,j) == max(fftUse(:,j)));
            peakVal(i,j) = fAx{i}(idxUse(peakIdx));
            
            % Create 3D matrix that holds different ranges
            rangeMatrix(i,j,k) = peakVal(i,j);

            % figure()
            % plot(fAx{i}(idxUse), fftUse(:,j)); 
            % set(gca,'YTickLabel',[]);
            % title('EEG Data: Song', i);
            % xlabel('Frequency (Hz)');
            % ylabel('Magnitude (A.U.)');
        end
    end
end


%% Get Difference in Hertz

% Actual Tempo (Hz) of Songs multiplied by 2, 4, and 8 due to how people
% can perceive a song as twice as fast or 4 times as fast.
original_tempo = [0.9328,1.1574,1.2376,1.3736,1.5244,1.6026,1.8116,2.0000,2.1368,2.5000];
tempo_mult2 = 2*[0.9328,1.1574,1.2376,1.3736,1.5244,1.6026,1.8116,2.0000,2.1368,2.5000];
tempo_mult4 = 4*[0.9328,1.1574,1.2376,1.3736,1.5244,1.6026,1.8116,2.0000,2.1368,2.5000];
tempo_mult8 = 8*[0.9328,1.1574,1.2376,1.3736,1.5244,1.6026,1.8116,2.0000,2.1368,2.5000];

% get the difference between participant max peak and target value
for k = 1:3
    for i = 1:10
        for j = 1:20
            tempo_peak2_diff(i,j,k) = rangeMatrix(i,j,k) - tempo_mult2(i);
            tempo_peak4_diff(i,j,k) = rangeMatrix(i,j,k) - tempo_mult4(i);
            tempo_peak8_diff(i,j,k) = rangeMatrix(i,j,k) - tempo_mult8(i);
        end
    end
end

%% Compare Matrices for Smallest Value

% Get smallest value out of the 3 peak_diff matrices
for k = 1:3
    for i = 1:10
        for j = 1:20
            % tempo_peak2_diff has value closest to 0
            if (abs(tempo_peak4_diff(i,j)) > abs(tempo_peak2_diff(i,j)))...
                    && (abs(tempo_peak2_diff(i,j))) < abs(tempo_peak8_diff(i,j))
                optimalPeak(i,j,k) = tempo_peak2_diff(i,j);
                predictedTempo(i,j,k) = rangeMatrix(i,j,k) / 2;

            % tempo_peak4_diff has value closest to 0
            elseif (abs(tempo_peak2_diff(i,j)) > abs(tempo_peak4_diff(i,j)))...
                    && (abs(tempo_peak4_diff(i,j))) < abs(tempo_peak8_diff(i,j))
                optimalPeak(i,j,k) = tempo_peak4_diff(i,j);
                predictedTempo(i,j,k) = rangeMatrix(i,j,k) / 4;

            % tempo_peak8_diff has value closest to 0
            elseif (abs(tempo_peak4_diff(i,j)) > abs(tempo_peak8_diff(i,j)))...
                    && (abs(tempo_peak8_diff(i,j))) < abs(tempo_peak2_diff(i,j))
                optimalPeak(i,j,k) = tempo_peak8_diff(i,j);
                predictedTempo(i,j,k) = rangeMatrix(i,j,k) / 8;
            end
        end
    end
end

%% Actual vs Predicted Tempo

for i = 1:10
    % This is a participant x frequency range matrix
    dataUse = squeeze(predictedTempo(i, :, :));
    diff_peak_tempo = dataUse - original_tempo(i);
    std(diff_peak_tempo);

    % Create boxplots showing predicted tempo (Hz) for 3 ranges
    figure();
    CategoricalScatterplot(dataUse, 'Labels',{'3-7','3-8','3-9'}...
        ,'BoxColor',[0.8471 0.8627 0.8392],'WhiskerColor',[0.8235 0.7412 0.0392]);
    title('Predicted Tempo in Range of Frequencies')
    xlabel('Range of Frequencies (Hz)')
    ylabel('Predicted Tempo (Hz)')
    legend([yline(original_tempo(i)),yline(2*original_tempo(i),'Color',[1,0,0])]...
        ,'Actual Tempo','Double Tempo',Location='northwest');
    ylim([0,5.0]);
end
%% Load Data

X = loadMultipleFiles();  % Feature vector

nComp = 1; % How many components to return
nReg = 7; % Regularization parameter - probably don't need to change ever 

%% Spatially Filter & Extract Useful Data

% Finds the most useful signal information in the dataset
% We only care about dataOut (data of the EEG that was filtered)
% Out of the 125 electric cusps used to capture brain waves, 
[dataOut,W,A,Rxx,Ryy,Rxy,dGen] = rcaRun125_parpoolAlready2021(...
    X, nReg, nComp, [], [], 0, []);

%% Convert to 2D Matrix & Plot

for i = 1:10
    % squeeze converts the 3D matrix into a 2D matrix
    % fft converts discrete signals from the time domain into a frequency domain
    dataFFT{i} = abs(fft(squeeze(dataOut{i})));
    dataMean{i} = abs(fft(mean(squeeze(dataOut{i}),2)));
    fAx{i} = computeFFTFrequencyAxis(size(dataFFT{i},1), 125);
    
    % Plot every participants' EEG data
    % figure();
    % plot(fAx{i}, dataFFT{i});
    % hold on    % Combines the plots
    % plot(fAx{i}, dataMean{i}, 'k', 'linewidth', 2);   % Add mean into the plot (black peaks)
    % xlim([0,15]);   % Limit the frequency from 0 to 15 Hz 
    % ylim([0,10*10^4]);
    % set(gca,'YTickLabel',[]);
    % title('EEG Data: Song', i);
    % xlabel('Frequency (Hz)');
    % ylabel('Magnitude (A.U.)');
end

%% Get Largest Peak Between 3 to (6+k) Hertz

% Using the largest peak values for ranges between 3-7 Hz, 3-8 Hz, & 3-9 Hz
% to analyze which range produces the best result

for k = 1
    for i = 9
        % finds index of values between 3 to (6+k) Hz
        idxUse = find(fAx{i} >= 0 & fAx{i} <= 15);
    
        % subset data values between 3 to (6+k) Hz
        fftUse = dataFFT{i}(idxUse,:);
    
        % create matrix of largest peak observed between 3 to (6+k) Hz
        for j = 1:20
            peakIdx = find(fftUse(:,j) == max(fftUse(:,j)));
            peakVal(i,j) = fAx{i}(idxUse(peakIdx));
            
            % Create 3D matrix that holds different ranges
            rangeMatrix(i,j,k) = peakVal(i,j);

            % figure()
            % plot(fAx{i}(idxUse), fftUse(:,j)); 
            % set(gca,'YTickLabel',[]);
            % title('EEG Data: Song', i);
            % xlabel('Frequency (Hz)');
            % ylabel('Magnitude (A.U.)');
        end
    end
end


%% Get Difference in Hertz

% Actual Tempo (Hz) of Songs multiplied by 2, 4, and 8 due to how people
% can perceive a song as twice as fast or 4 times as fast.
original_tempo = [0.9328,1.1574,1.2376,1.3736,1.5244,1.6026,1.8116,2.0000,2.1368,2.5000];
tempo_mult2 = 2*[0.9328,1.1574,1.2376,1.3736,1.5244,1.6026,1.8116,2.0000,2.1368,2.5000];
tempo_mult4 = 4*[0.9328,1.1574,1.2376,1.3736,1.5244,1.6026,1.8116,2.0000,2.1368,2.5000];
tempo_mult8 = 8*[0.9328,1.1574,1.2376,1.3736,1.5244,1.6026,1.8116,2.0000,2.1368,2.5000];

% get the difference between participant max peak and target value
for k = 1:3
    for i = 1:10
        for j = 1:20
            tempo_peak2_diff(i,j,k) = rangeMatrix(i,j,k) - tempo_mult2(i);
            tempo_peak4_diff(i,j,k) = rangeMatrix(i,j,k) - tempo_mult4(i);
            tempo_peak8_diff(i,j,k) = rangeMatrix(i,j,k) - tempo_mult8(i);
        end
    end
end

%% Compare Matrices for Smallest Value

% Get smallest value out of the 3 peak_diff matrices
for k = 1:3
    for i = 1:10
        for j = 1:20
            % tempo_peak2_diff has value closest to 0
            if (abs(tempo_peak4_diff(i,j)) > abs(tempo_peak2_diff(i,j)))...
                    && (abs(tempo_peak2_diff(i,j))) < abs(tempo_peak8_diff(i,j))
                optimalPeak(i,j,k) = tempo_peak2_diff(i,j);
                predictedTempo(i,j,k) = rangeMatrix(i,j,k) / 2;

            % tempo_peak4_diff has value closest to 0
            elseif (abs(tempo_peak2_diff(i,j)) > abs(tempo_peak4_diff(i,j)))...
                    && (abs(tempo_peak4_diff(i,j))) < abs(tempo_peak8_diff(i,j))
                optimalPeak(i,j,k) = tempo_peak4_diff(i,j);
                predictedTempo(i,j,k) = rangeMatrix(i,j,k) / 4;

            % tempo_peak8_diff has value closest to 0
            elseif (abs(tempo_peak4_diff(i,j)) > abs(tempo_peak8_diff(i,j)))...
                    && (abs(tempo_peak8_diff(i,j))) < abs(tempo_peak2_diff(i,j))
                optimalPeak(i,j,k) = tempo_peak8_diff(i,j);
                predictedTempo(i,j,k) = rangeMatrix(i,j,k) / 8;
            end
        end
    end
end

%% Actual vs Predicted Tempo

for i = 1:10
    % This is a participant x frequency range matrix
    dataUse = squeeze(predictedTempo(i, :, :));
    diff_peak_tempo = dataUse - original_tempo(i);
    std(diff_peak_tempo);

    % Create boxplots showing predicted tempo (Hz) for 3 ranges
    figure();
    CategoricalScatterplot(dataUse, 'Labels',{'3-7','3-8','3-9'}...
        ,'BoxColor',[0.8471 0.8627 0.8392],'WhiskerColor',[0.8235 0.7412 0.0392]);
    title('Predicted Tempo in Range of Frequencies')
    xlabel('Range of Frequencies (Hz)')
    ylabel('Predicted Tempo (Hz)')
    legend([yline(original_tempo(i)),yline(2*original_tempo(i),'Color',[1,0,0])]...
        ,'Actual Tempo','Double Tempo',Location='northwest');
    ylim([0,5.0]);
end
