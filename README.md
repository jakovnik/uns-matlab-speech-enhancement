# uns-matlab-speech-enhancement
Matlab library for speech enhancement
## Structure of the repository
- noise_est_min_stat.m - It estimates noise power spectrum density using algorithm proposed in Martin, R. "Noise power spectral density estimation based on optimal smoothing and minimum statistics" _IEEE Transactions on Speech and Audio Processing_, 2001, 9, 504-512.
- noise_est_snr_recursive.m - It estimates noise power spectrum density using algorithm proposed in Lin, L.; Holmes, W. H. & Ambikairajah, E. "Adaptive noise estimation algorithm for speech enhancement" _Electron. Lett._, 2003, 39, 754-755.
- load_audio_from_dir.m - It loads audio files from test and referent folder and store them into a structure which will be used in testing steps.
- signal_segmentation.m - It segments signal into frames.
- evaluate_noise_estimation_error.m - It evaluates quality of noise estimation by comparison of estimated and referent noise level.
- test_noise_min_stat.m .m - It is a scrpit for evaluation of noise_est_min_stat on NOIZEUS database. It is also an example how to use function noise_est_min_stat.m 
- test_noise_est_snr_recursive.m - It is a scrpit for evaluation of noise_est_snr_recursive on NOIZEUS database. It is also an example how to use function noise_est_snr_recursive.m 
 



