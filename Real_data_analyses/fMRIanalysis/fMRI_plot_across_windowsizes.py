import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io as sio

# Load data
data = np.load('pid_vals_allwindows_shrink_random_shortdelay.npz')
syns = data['Nsyns_measures_tot']  # shape: [triplets, subjects, orders, methods, time_windows]
reds = data['Nreds_measures_tot']
uniques = data['Nuniques_measures_tot']
joints = data['Njoints_measures_tot']
time_windows = data['time_windows']
measures = ['Nsyns_measures_tot', 'Nreds_measures_tot', 'Nuniques_measures_tot', 'Njoints_measures_tot']
clean_matrix = np.zeros((15, 100, 20, 4, 4))

print(syns.shape)
for idx, key in enumerate(measures):
    clean_matrix[:15, :, :, :, idx] = data[key][:15, :, :20, :, -1]
sio.savemat("fMRI_alltimewindows.mat", {
    'results': clean_matrix,
})

# Configuration
measures = {'Joint': joints, 'Synergy': syns, 'Redundancy': reds, 'Unique': uniques}
methods = ['NoBias', 'InfoBias', 'ShuffSub', 'Venkatesh']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
orders = np.arange(1, syns.shape[2]+1)  # x-axis: dimensions
custom_xticks = [1, 2] + list(range(4, 21, 2))
# Create output directory
os.makedirs("pid_plots", exist_ok=True)

# Plot for each time window
for t_idx, t_win in enumerate(time_windows):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for i, (measure_name, data_array) in enumerate(measures.items()):
        ax = axes[i]

        for method_idx in range(data_array.shape[3]):
            # Extract data: [triplets, subjects, orders]
            vals = data_array[:, :, :, method_idx, t_idx]

            # Mean and SEM across triplets and subjects
            mean_vals = np.nanmean(vals, axis=(0, 1))
            sem_vals = np.nanstd(vals, axis=(0, 1)) / np.sqrt(np.sum(~np.isnan(vals), axis=(0, 1)))
            error = 2 * sem_vals  # 2 * SEM

            ax.plot(orders[:20], mean_vals[:20], label=methods[method_idx], color=colors[method_idx])
            ax.fill_between(orders[:20], mean_vals[:20] - error[:20], mean_vals[:20] + error[:20], alpha=0.3, color=colors[method_idx])

        ax.set_title(measure_name)
        ax.set_xlabel('Dimension d')
        # ax.set_xlim([1, 20])
        ax.set_ylabel('Value')
        ax.set_xticks(custom_xticks)
        ax.grid(True)
        ax.legend()

    plt.tight_layout()
    plt.suptitle(f'PID Components - Time Window: {t_win}', fontsize=16, y=1.02)

    # Save figure
    fname = f'pid_plots/pid_plot_timewin_{t_win}.png'
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fname}")
