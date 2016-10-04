function y = signal_segmentation(sig, window_fnc, window_shift)
% Split signal into frames.
% The first frame is centered at window_shift
%
% Inputs:
%   sig - input signal
%   window_fnc - window function
%   window_shift - window shift [in samples]
%
% Output:
%   y - segmented signal - each column is a frame
%
% Author: jakovnik@gmail.com
% Date: 2016/8/31
%

sig = sig(:);
window_fnc = window_fnc(:);
wds = length(window_fnc);
nSegments = floor(length(sig)/window_shift) + 1;
sig(nSegments*window_shift+ceil(0.5*wds)) = 0;
sig = [zeros(ceil(0.5*wds-window_shift),1);sig];
y = sig([1:wds]'* ones(1, nSegments) + ...
  ones(wds,1) * [0:nSegments-1] * window_shift);
clear sig
for frm_idx = 1:nSegments
   y(:,frm_idx) = window_fnc.*y(:,frm_idx);
end
