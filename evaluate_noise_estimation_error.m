function [mse, med_se, d] = evaluate_noise_estimation_error(estimated_psd, true_psd, ignore_n_last_frames)
% It calculates objective measures for evaluatin accuracy of the
% noise-estimation algorithms, i.e. mean squared error (mse) and
% median squared error (med_se)
%
% Inputs:
%  estimated_psd - estimated noise power spectral density
%  true_psd - true noise power spectral density
%  ignore_n_last_frames - number of the last framest to ignore. Default
%       value is zero
% In these matrices the first index corresponds to the frequency, and
% the second index to the time. The input matrices shoud have the same
% dimensions.
%
% Outputs:
%  mse - mean square error
%  med_se - median square error
%  d - psd square error over time
%
% Author: jakovnik@gmail.com
% Date: 2016/8/26
%

if nargin < 3
   ignore_n_last_frames = 0;
end

%% Check dimensions match
[n1, m1] = size(estimated_psd);
[n2, m2] = size(true_psd);
if or(n1 ~= n2, m1 ~= m2)
  error('Input data dimensions mismatch.')
end
estimated_psd = estimated_psd(:,1:end-ignore_n_last_frames);
true_psd = true_psd(:,1:end-ignore_n_last_frames);
%% Evaluate objective measures
d = sum((estimated_psd - true_psd).^2)./sum(true_psd.^2);
mse = mean(d);
med_se = median(d);
