%% Set parameters
%! Set path to the NOIZEUS database (do not forget '/' at end of the path)
%database_path = 'c:/Users/Nikša/Documents/doas/speech_enhancement/noizeus/';
%database_path = '/mnt/users/jakovljevic/niksa/_working/speech_enhancement/noizeus/'
%! Set the name of the output .csv file - (it writes in append mode)
% output_file = 'noise_est_min_stats.csv';
sample_rate = 8000;
window_duration = 32e-3;
window_shift = 8e-3;
param.D = 96;
param.V = 12;
param.tshift = window_shift;
param.tdecay = 64e-3;
param.alpha_max = 0.96;
param.alpha_min = 0.30;
param.alpha_c_min = 0.7;
param.min_y_energy = 1e-9;
param.inv_q_max = 1/2;
param.beta_max = 0.8;
%% Create the list of files for testing on NOIZEUS database
test_set = cell(32,1);
test_set{1}.row = 2; test_set{1}.type = 'airport'; test_set{1}.snr = '15';
test_set{2}.row = 3; test_set{2}.type = 'airport'; test_set{2}.snr = '10';
test_set{3}.row = 4; test_set{3}.type = 'airport'; test_set{3}.snr = '5';
test_set{4}.row = 5; test_set{4}.type = 'airport'; test_set{4}.snr = '0';
test_set{5}.row = 6; test_set{5}.type = 'babble'; test_set{5}.snr = '15';
test_set{6}.row = 7; test_set{6}.type = 'babble'; test_set{6}.snr = '10';
test_set{7}.row = 8; test_set{7}.type = 'babble'; test_set{7}.snr = '5';
test_set{8}.row = 9; test_set{8}.type = 'babble'; test_set{8}.snr = '0';
test_set{9}.row = 10; test_set{9}.type = 'car'; test_set{9}.snr = '15';
test_set{10}.row = 11; test_set{10}.type = 'car'; test_set{10}.snr = '10';
test_set{11}.row = 12; test_set{11}.type = 'car'; test_set{11}.snr = '5';
test_set{12}.row = 13; test_set{12}.type = 'car'; test_set{12}.snr = '0';
test_set{13}.row = 14; test_set{13}.type = 'exhibition'; test_set{13}.snr = '15';
test_set{14}.row = 15; test_set{14}.type = 'exhibition'; test_set{14}.snr = '10';
test_set{15}.row = 16; test_set{15}.type = 'exhibition'; test_set{15}.snr = '5';
test_set{16}.row = 17; test_set{16}.type = 'exhibition'; test_set{16}.snr = '0';
test_set{17}.row = 18; test_set{17}.type = 'restaurant'; test_set{17}.snr = '15';
test_set{18}.row = 19; test_set{18}.type = 'restaurant'; test_set{18}.snr = '10';
test_set{19}.row = 20; test_set{19}.type = 'restaurant'; test_set{19}.snr = '5';
test_set{20}.row = 21; test_set{20}.type = 'restaurant'; test_set{20}.snr = '0';
test_set{21}.row = 22; test_set{21}.type = 'station'; test_set{21}.snr = '15';
test_set{22}.row = 23; test_set{22}.type = 'station'; test_set{22}.snr = '10';
test_set{23}.row = 24; test_set{23}.type = 'station'; test_set{23}.snr = '5';
test_set{24}.row = 25; test_set{24}.type = 'station'; test_set{24}.snr = '0';
test_set{25}.row = 26; test_set{25}.type = 'street'; test_set{25}.snr = '15';
test_set{26}.row = 27; test_set{26}.type = 'street'; test_set{26}.snr = '10';
test_set{27}.row = 28; test_set{27}.type = 'street'; test_set{27}.snr = '5';
test_set{28}.row = 29; test_set{28}.type = 'street'; test_set{28}.snr = '0';
test_set{29}.row = 30; test_set{29}.type = 'train'; test_set{29}.snr = '15';
test_set{30}.row = 31; test_set{30}.type = 'train'; test_set{30}.snr = '10';
test_set{31}.row = 32; test_set{31}.type = 'train'; test_set{31}.snr = '5';
test_set{32}.row = 33; test_set{32}.type = 'train'; test_set{32}.snr = '0';
% Prepare file for writing
fex = exist(output_file,'file');
fileID = fopen(output_file,'at');
if fex == 0 %print only first time
  fprintf(fileID,['%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'...
    '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'],...
    'type', 'snr', 'window_duration', 'window_shift',...
    'D','V','tdecay','alpha_max','alpha_min',...
    'alpha_c_min','min_y_energy','inv_q_max','beta_max',...
    'mse', 'med se', 'std', 'iqr','prc5','prc25','prc75','prc95','max',...
    'comment', 'max mse');
end
win_dur_smp = round(window_duration*sample_rate);
window_fnc = hamming(win_dur_smp);
L = 2^ceil(log2(win_dur_smp));
win_sft_smp = round(window_shift*sample_rate);
for j = 1:length(test_set)
  type = test_set{j}.type;
  snr = test_set{j}.snr;
  test_dir = [database_path,type,'/',snr,'dB'];
  ref_dir = [database_path,'clean'];
  %% Initialize
  test_files = load_audio_from_dir(test_dir,ref_dir);
  n_files = length(test_files);
  %% Process files
  all_d = [];
  the_worst = '';
  max_mse = -1;
  for k = 1:n_files
    sig = test_files{k}.sig;
    s_sig = signal_segmentation(sig,window_fnc,win_sft_smp);
    s_noise = signal_segmentation(test_files{k}.noise,window_fnc,win_sft_smp);
    sig_psd = fft(s_sig,L);
    sig_psd = abs(sig_psd(1:L/2+1,:)).^2;
    noise_psd = fft(s_noise,L);
    noise_psd = abs(noise_psd(1:L/2+1,:)).^2;
    [sigma_N_2, sn] = noise_est_min_stat(sig_psd,param);
    [mse, med_se, d] = evaluate_noise_estimation_error(sigma_N_2, noise_psd,2);
    disp([test_files{k}.name,' ',num2str(mse),' ',num2str(med_se)])
    all_d = [all_d, d];
    if mse > max_mse
      max_mse = mse;
      the_worst = test_files{k}.name;
    end
  end
  prc = prctile(all_d, [5; 25; 50; 75; 95]);
  mse = mean(all_d)
  med_se = prc(3)
  std_ = std(all_d)
  iqr_ = prc(4)-prc(2);
  fprintf(fileID,['%s\t%s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t'...
    '%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%s\t%e\n'],...
    type, snr, window_duration, window_shift,...
    param.D,param.V,param.tdecay,param.alpha_max,...
    param.alpha_min,param.alpha_c_min,param.min_y_energy,...
    param.inv_q_max, param.beta_max, mse, med_se, std_, iqr_,...
    prc(1),prc(2),prc(4),prc(5),max(all_d),...
    the_worst, max_mse);  
end
fclose(fileID);

