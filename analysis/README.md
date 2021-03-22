# respirationCA - Data analysis

## Software requirements

- R
- Matlab
- Kubios
- Helper functions<sup>1</sup>

<sup>1</sup>Part of repository (see "assets/")

## Behavioral data analysis

Directory: behavior/

1. respirationCA_trials.txt - behavioral data (1 trial each row)
2. behav.R - main script to preprocess behavioral data
3. plot_detection.R - plot detection rates
4. plot_conf.R - plot confidence ratings
5. plot_resp_time.R - plot detection and confidence response times
6. plot_resp_t_conf.R - plot detection and confidence response times split by confidence

## ECG data analysis

Directory: ecg/

1. preprocessing with Kubios to get R peaks
2. ecg_stim2peak.R - stimulus onset relative to previous R peak and cardiac cycle -> ecgdata.csv
3. test_uniform_stim.R - Rayleigh tests for uniform distribution of stimulus onsets -> ecg_final.csv
4. bin_resp1.R - test if hit rate differs between 4 time bins 0.0-0.8 s after R peak
5. plot_mean_angle.R - plot distribution of mean angles for hits, misses and correct rejections

## Respiration data analysis

Directory: respiration/

1. resp_analysis.m - main script to load/preprocess data, find cycles, and locate stimulus onsets
2. merge_behav_resp.m - merges output of resp_analysis.m with behavior data
3. circular_distribution.R - test and plot circular distribution of stimulus-response conditions
4. entrainment.R - checks circular distance to block mean angle
5. breath_duration.R - analyze and plot respiration cycle durations
6. plv.R - calculate phase-locking values for cardiac and respiratory angles

## Oximeter data analysis

Directory: oxi/

1. oxi_analysis.m - main script to load, preprocess, analyze, and plot data

## Somatosensory evoked potential (SEP) analysis

Directory: respiration/

1. sep_loop.m - main script to run sep_analysis.m for all participants and plot results
