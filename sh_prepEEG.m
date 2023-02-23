function EEG = sh_prepEEG(Data, filename, chanlocs, fs)
    EEG = eeg_emptyset();
    EEG.setname = filename;
    EEG.filename = filename;
    EEG.nbchan = size(Data, 1);
    EEG.trials = 1;
    EEG.pnts = size(Data, 2);
    EEG.srate = fs;
    EEG.xmin = 0;
    EEG.xmax = EEG.pnts/EEG.srate;
    EEG.times = 1/EEG.srate:1/EEG.srate:EEG.xmax;
    EEG.data = Data;
    EEG.chanlocs = chanlocs;
end