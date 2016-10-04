function test_files = load_audio_from_dir(test_dir,ref_dir)
% It load audio files from test and referent folder and store them into a
% structure. Each element of the structure contains the name of test file,
% samples of test file and noise samples (test-referent samples).
%
% Inputs:
%   test_dir - path to the folder containing test files
%   ref_dir - path to the folder containing corresponding clean audio files
%
% Output:
%   test_files - structure containing audio samples of noisy signal (sig),
%     noise (noise), and file name (name).
%
% Author: jakovnik@gmail.com
% Date: 2016/8/26

%% Load test files (noisy speech)
wav_list = dir([test_dir,'/*.wav']);
n_files = length(wav_list);
ts_list = cell(n_files,1);
for k = 1:n_files;
  fname = wav_list(k).name;
  vname = matlab.lang.makeValidName(fname(1:end-4));
  ts_list{k} = vname;
  command = strcat(vname, ' = audioread( '' ', test_dir, '/' ,fname, ''');');
  eval( command );
end
%% Load referent files (clean speech)
wav_list = dir([ref_dir,'/*.wav']);
for k = 1:length(wav_list);
  fname = wav_list(k).name;
  vname = matlab.lang.makeValidName(fname(1:end-4));
  command = strcat(vname, ' = audioread( '' ', ref_dir, '/' ,fname, ''');');
  eval( command );
end
%% Store data
test_files = cell(n_files,1);
for k = 1:n_files
  if exist(ts_list{k}(1:4),'var') ~= 1
    warning('load_audio_from_dir:WAR1',...
      '%s does not exist - skip file %s', ts_list{k}(1:4),ts_list{k});
  end
  b.name = ts_list{k};
  eval(['b.sig = ',ts_list{k},';']);
  eval(['b.noise = ',ts_list{k},' - ',ts_list{k}(1:4),';']);
  test_files{k} = b;
end
