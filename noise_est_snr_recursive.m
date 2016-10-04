function [sigma_N_2_final, sigma_N_2, alpha_trace] = noise_est_snr_recursive(abs_Y_2, param)
% Estimate noise PSD using algorithm proposed in [1]

% Inputs:
%   abs_Y_2 - power spectrum distribution of noisy signal
%             Instead of PSD, any filter bank output can be used
%   param - parameters of the algorithm
%       + window_length - the length of the analysis window in frames
%       + use_median - a switch - whether to use median or average 
%             to estimate noise PSD
%       + a - parameter in sigmoid function which defines slope
%       + T - parameter in sigmodi function which defines symmetry
%          alpha(snr) = 1 / (1 + exp(-a (snr-T)))
%       + alpha_final - one pole psd smoothing filter 

% [1] Lin, L.; Holmes, W. H. & Ambikairajah, E. 
%    "Adaptive noise estimation algorithm for speech enhancement"
%    Electron. Lett., 2003, 39, 754-755

% author:  jakovnik@gmail.com
% date: 2016/08/17

if or(nargin < 1, nargin > 2)
  error ('invalid number of input arguments');
end
%% Set default values
if nargin == 1
  param = [];
end
L = 10; %value between 5 and 10.
if isfield(param,'window_length')
  L = param.window_length;
end
use_med = 0; 
if isfield(param,'use_median')
  use_med = param.use_median;
end
a = 20; %value between 15 and 30
if isfield(param,'a')
  a = param.a;
end
T = 1.5; % value around 1.5
if isfield(param,'T')
  T = param.T;
end
alpha_fin = 0.8; %value between 0.5 and 0.95
if isfield(param,'alpha_final')
  alpha_fin = param.alpha_final;
end

%% Initialize
[n_f_bins, n_frames] = size(abs_Y_2);
sigma_N_2 = zeros(n_f_bins, n_frames);
sigma_N_2_final = zeros(n_f_bins, n_frames);
alpha_trace = zeros(n_f_bins,n_frames);
sigma_N_2(:,1) = abs_Y_2(:,1);
sigma_N_2_final(:,1) = sigma_N_2(:,1);

%% Process
for fidx = 2:n_frames
  if (use_med)                                         % par. after eq. (4)
    mean_sigma_N_2 = median(sigma_N_2(:,max(fidx-L,1):fidx-1),2);
  else
    mean_sigma_N_2 = mean(sigma_N_2(:,max(fidx-L,1):fidx-1),2);
  end
  snr = abs_Y_2(:,fidx)./mean_sigma_N_2;               % par. after eq. (4)
  alpha = 1./(1 + exp(-a*(snr-T)));                               % eq. (4)
  sigma_N_2(:,fidx) = alpha.*sigma_N_2(:,fidx-1) + ...
                      (1 - alpha).*abs_Y_2(:,fidx);               % eq. (3)
  sigma_N_2_final(:,fidx) = alpha_fin*sigma_N_2_final(:,fidx-1) +...
                            (1 - alpha_fin).*sigma_N_2(:,fidx);   % eq. (5)
  alpha_trace(:,fidx) = alpha;
end