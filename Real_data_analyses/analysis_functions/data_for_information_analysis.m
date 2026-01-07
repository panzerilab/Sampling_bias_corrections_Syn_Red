function data_for_information_analysis
%%
    original_data                   = dir([pwd,'/Data/*.mat']);
    mkdir('Data/data_for_information_analysis/');
    for ind                        =  1:length(original_data)
        ind
        data_mouse              = load(['Data/',original_data(ind).name]); %%% Load all the data
        Behavior                = data_mouse.Behavior;                        %%% Load only trials during the task
        % Un-filtered DFF
        DFF_corrected_Beh       = data_mouse.Behavior.Fluorescence.DeltaFCellTrialsCorrected;   %%% We consider a subset of task trials 
        % Sgolay filter DFF
        DFF_sgolay_Beh          = smoothdata(DFF_corrected_Beh, 1, 'sgolay', [10 10], 'Degree', 5);
        % Butterworth filter DFF
        DFF_smoothed_Beh        = data_mouse.Behavior.Fluorescence.DeltaFCellTrialsCorrectedSmooth;   %%% We consider a subset of task trials
        % FIR filter DFF
        DFF_FIR_Beh             = zeros(size(data_mouse.Behavior.Fluorescence.DeltaFCellTrialsCorrected,1),size(data_mouse.Behavior.Fluorescence.DeltaFCellTrialsCorrected,2),size(data_mouse.Behavior.Fluorescence.DeltaFCellTrialsCorrected,3));
        SGPoly                  = polyfit(-10:10, [zeros(1,10),1, zeros(1,10)],5);
        FIR_Filt                = fliplr(polyval(SGPoly,-10:10));
        for trial               = 1:size(data_mouse.Behavior.Fluorescence.DeltaFCellTrialsCorrected,2)
            for neurons         = 1:size(data_mouse.Behavior.Fluorescence.DeltaFCellTrialsCorrected,3)
                DFF_FIR_Beh(:,trial,neurons)              = conv(data_mouse.Behavior.Fluorescence.DeltaFCellTrialsCorrected(:,trial,neurons), FIR_Filt, 'same');  
            end
        end
        responsivecells         = Behavior.Fluorescence.responsiveCells==1;   %%% Extract responsive cells
        brightcells             = Behavior.Fluorescence.brightcells'==1;   %%% Extract brightcells
        n_responsive_neurons    = size(DFF_corrected_Beh(:,:,responsivecells),3);    %%% Responsive neurons
        timepoints              = size(DFF_corrected_Beh(:,:,responsivecells),1);    %%% Size of timepoints
        name_exp_mice           = Behavior.Animal;                            %%% Name of the individual
        date_exp_mice           = Behavior.Date;                              %%% Day of the experiment
        cell_location_mice      = Behavior.Fluorescence.CellLoc;              %%% XY coordinates of all neurons per experiment
        stimonset               = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'StimOnset')))); %%% A vector of time frames when a stimulus occured in different trials (at tf=31)
        StimFrequency           = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'StimFrequency')))); %%% A vector of frequencies in different trials (at tf=31)
        early_trials            = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'Early')))); %%% A vector of early trials (at tf=31, 1==early trials, 0==otherwise)
        early_FA_trials         = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'EarlyFalseAlarm')))); %%% A vector of early false alarms (at tf=31, 1==early false alarms, 0==otherwise)
        early_hit_trials        = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'EarlyHit')))); %%% A vector of early hits (at tf=31, 1==early hits, 0==otherwise)
        hit_trials              = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'Hit')))); %%% A vector of hit trials (at tf=31, 1==hit trials, 0==otherwise)
        miss_trials             = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'Miss')))); %%% A vector of miss trials (at tf=31, 1==miss trials, 0==otherwise)
        FA_trials               = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'FalseAlarm')))); %%% A vector of false alarms (at tf=31, 1==false alarms, 0==otherwise)
        CR_trials               = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'CorrectReject')))); %%% A vector of correct rejection trials (at tf=31, 1==correct rejection trials, 0==otherwise)
        FirstResponse           = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'FirstResponse')))');        %find where the first lick occured relative to the image frame timing
        save(['Data/data_for_information_analysis/processed_for_information_analysis_',char(original_data(ind).name)],...
            'DFF_corrected_Beh','DFF_smoothed_Beh','DFF_FIR_Beh',...
            'DFF_sgolay_Beh','stimonset','StimFrequency','early_trials','early_FA_trials', 'early_hit_trials',...
            'hit_trials','miss_trials','FA_trials','CR_trials','responsivecells','brightcells', 'name_exp_mice','date_exp_mice',...
            'cell_location_mice','n_responsive_neurons','timepoints','FirstResponse');    %%% Save (what I need) to Data/data_for_information_analysis folder
    end
end