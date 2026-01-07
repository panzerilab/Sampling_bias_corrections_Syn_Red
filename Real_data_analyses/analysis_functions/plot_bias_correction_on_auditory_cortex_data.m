function plot_bias_correction_on_auditory_cortex_data(Rbins, type_data, deconv_type, bias_correction_ALL)
%%
    bin_size = 10; type_trials = 'ALL';    
    load(['Results/Real_data_analysis/info_results_pairwise_PEAK_PID_sliding_window_',num2str(bin_size),'tp_ALLtrials_',char(bias_correction_ALL),'_',num2str(Rbins),'.mat'],'info_results_pairwise_PEAK_PID');
    plot_ALL_correct_incorrect_auditory_cortex_data(info_results_pairwise_PEAK_PID,type_trials, bias_correction_ALL, bin_size,type_data, deconv_type);
    
end