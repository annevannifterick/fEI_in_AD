% To calculate spectrum DFA and fEI for all recordings

% Colors (1st row: TDC; 2nd row: STXBP1)
dark = [0.3255, 0.3176, 0.3294; 0.2235, 0.4157, 0.6941];
light = [0.5020, 0.5216, 0.5216; 0.4471, 0.5765, 0.7961];

% Paths
path_input = 'data_preproc';
path_output = 'EIspec';

if ~exist(path_output, 'dir')
    mkdir(path_output);
end

addpath("support_func")

% Specify parameters

% Parameters for frequencies
freqRes = 1;    % Hz
freqMin = 1;    % Hz
freqMax = 45;   % Hz
freqs = freqMin:freqRes:freqMax-1;
nFreqs = length(freqs);

% Parameters for DFA and fE/I
overlap_DFA = 0.5;
fitRange_DFA = [4, 20];
overlap_EI = 0.8;

% Index recordings in the input folder
recordings = dir(fullfile(path_input, 'case*.mat'));

% Compute DFA and fE/I spectrum per recording
for recordingIdx = 1 : length(recordings)
    disp(['Recording: ', num2str(recordingIdx)]);
    
    chanlocs = [];  % Channel locations if  you want to plot the results on a topoplot
    
    % Load the data
    load(fullfile(recordings(recordingIdx).folder, recordings(recordingIdx).name));
    
    % Sampling rate (modify if it differs per recording!)
    fs = 1250;  % Sampling rate in Hz
    windowSize_EI = 5 * fs;
    
    % Transpose signal
    data = signal_cut';
    
    % Number of channels
    nChannels = size(data, 1);
    
    % Prepare matrices for DFA and fE/I
    DFA_spectrum = nan(nChannels, nFreqs);
    EI_spectrum = nan(nChannels, nFreqs);
    
    for fIdx = 1 : nFreqs
        disp(['Frequency bin: ', num2str(freqs(fIdx)), '-', num2str(freqs(fIdx)+freqRes), ' Hz']);
        
        % FIR filter from EEGLAB
        EEG = sh_prepEEG(data, '', chanlocs, fs);
        EEG_filt = pop_eegfiltnew(EEG, freqs(fIdx), freqs(fIdx)+freqRes);
        
        % Get amplitude envelope using the Hilbert transform
        ampEnv = abs(hilbert(EEG_filt.data'));
        
        % Compute DFA for frequency bin
        DFA_spectrum(:, fIdx) = calculateDFA(ampEnv, fitRange_DFA, [0.5, 100], EEG_filt.srate, overlap_DFA)';
        
        % Compute fE/I for frequency bin
        %[~, EI_spectrum(:, fIdx)] = calculateFEI(ampEnv, windowSize_EI, overlap_EI, DFA_spectrum(:, fIdx));
        [EI_spectrum(:, fIdx), ~] = calculateFEI(ampEnv, windowSize_EI, overlap_EI, DFA_spectrum(:, fIdx));
    end
    
    % Save 
    filename_output = strrep(fullfile(recordings(recordingIdx).folder, recordings(recordingIdx).name), path_input, path_output);
    save(filename_output, 'DFA_spectrum', 'EI_spectrum', 'freqRes', 'freqMin', 'freqMax', 'freqs', 'overlap_DFA', 'fitRange_DFA', 'overlap_EI', 'windowSize_EI', 'nChannels');
end