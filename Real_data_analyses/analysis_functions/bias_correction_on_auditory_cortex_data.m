function bias_correction_on_auditory_cortex_data(Rbins, spike_train_alltrials, deconv_type, type_data, bias_correction_ALL,redundancy_measure)
    %%
    info_type = 'stimulus_information'; type_trials = 'ALL';
    [info_results_pairwise_PEAK_PID, bin_size] = pairwise_neuron_PEAK_information_PID(Rbins, info_type,type_trials,spike_train_alltrials, deconv_type, type_data, bias_correction_ALL, redundancy_measure);
    save(['Results/Real_data_analysis/info_results_pairwise_PEAK_PID_sliding_window_',num2str(bin_size),'tp_ALLtrials_',char(bias_correction_ALL),'_',num2str(Rbins),'.mat'],'info_results_pairwise_PEAK_PID');

end