function [sigma_N_2, P, alpha] = noise_est_min_stat(abs_Y_2, param)
% Estimate noise spectrum using minimum statistics as proposed in [1]
%
% Inputs:
%   abs_Y_2 - noisy signal PSD (one column per frame)
%   param - parameters of the algorithm
%       + D - number of frames in analysing window for minimum search
%       + V - number of frames in analysing subwindow 
%       + tshift - frame shift in seconds
%       + tdecay - decay time from signal peak to the noise level
%       + alpha_max - maximum value of smoothing factor
%       + alpha_min - minimum value of smoothing factor 
%       + alpha_c_min - minimum value of alpha_c
%       + min_y_energy - minimum energy
%       + inv_q_eq_max - maximum value of Q_eq^-1
%       + beta_max / maximum value of beta
%       + noise_slope_max_table - for comparison of Q^-1
%       + a_v - factor used in B_c calcluation
%
% Outputs:
%  sigma_N_2 - estimated noise PSD
%  P - smoothed periodogram
%  alpha - optimal smoothing factor

% [1] Martin, R. 
%     "Noise power spectral density estimation based on optimal smoothing
%      and minimum statistics" 
%      IEEE Transactions on Speech and Audio Processing, 2001, 9, 504-512

% author:  jakovnik@gmail.com
% date: 2016/08/23

if or(nargin < 1, nargin > 2)
  error ('invalid number of input arguments');
end
%% Set default values
if nargin == 1
  param = [];
end
D = 96;
if isfield(param,'D')
  D = param.D;
end
V = 12; 
if isfield(param,'V')
  V = param.V;
end
tshift = 16e-3; 
if isfield(param,'tshift')
  tshift = param.tshift;
end
tdecay = 64e-3;
if isfield(param,'tdecay')
  tdecay = param.tdecay;
end
alpha_max = 0.96;
if isfield(param,'alpha_max')
  alpha_max = param.alpha_max;
end
alpha_min = 0.30;
if isfield(param,'alpha_min')
  alpha_min = param.alpha_min;
end
alpha_c_min = 0.7;
if isfield(param,'alpha_c_min')
  alpha_c_min = param.alpha_c_min;
end
min_y_energy = 1e-9;
if isfield(param,'min_y_energy')
  min_y_energy = param.min_y_energy;
end
inv_q_eq_max = 1/2;
if isfield(param,'inv_q_eq_max')
  inv_q_eq_max = param.inv_q_eq_max;
end
beta_max = 0.8;
if isfield(param,'beta_max')
  beta_max = param.beta_max;
end
noise_slope_max_table = [0.03 8; 
                         0.05 4; 
                         0.06 2];
if isfield(param,'noise_slope_max_table')
  noise_slope_max_table = param.noise_slope_max_table;
end
a_v = 2.12;
if isfield(param,'a_v')
  a_v = param.a_v;
end
MAX_VALUE = 1e9;

U = round(D/V);
M_D = mean_of_the_min_table3(D);
M_V = mean_of_the_min_table3(V);
[L,n_frames] = size(abs_Y_2);
snr_exp = -tshift/tdecay;

% initialize
sigma_N_2 = zeros(L,n_frames);
P = zeros(L,n_frames);
alpha = zeros(L,n_frames);
% calculation
buff_idx = 1;
subwc = V;
P_prev = abs_Y_2(:,1);
alpha_c_prev = 1; %set initial value (it is suppose noise at begin)
sigma_N_2_prev = abs_Y_2(:,1);
sigma_N_2(:,1) = abs_Y_2(:,1);
mean_P = abs_Y_2(:,1);
mean_P_2 = mean_P.^2;
actmin_buff = MAX_VALUE*ones(L,U);
actmin = MAX_VALUE*ones(L,1);
actmin_sub = MAX_VALUE*ones(L,1);
lmin_flag = zeros(L,1);
for lambda = 1:n_frames
  % compute the smoothing parameter alpha
  y_energy = sum(abs_Y_2(:,lambda));
  p_energy = sum(P_prev);
  if y_energy == 0
    y_energy = min_y_energy; 
  end
  alpha_c = 1/(1 + (p_energy/y_energy - 1)^2);                      %eq.  9
  alpha_c = 0.7*alpha_c_prev + 0.3*max(alpha_c,alpha_c_min);        %eq. 10
  alpha_c_prev = alpha_c;
  alpha_opt = alpha_max * alpha_c ./...
    (1 + (P_prev./sigma_N_2_prev - 1).^2);                          %eq. 11
  snr = p_energy/sum(sigma_N_2_prev);
  alpha_opt = max(alpha_opt,min(alpha_min,snr^snr_exp));            %eq.+12
  alpha(:,lambda) = alpha_opt; 
  % compute smoothed power P(\lambda,k)
  P(:,lambda) = alpha_opt .* P_prev + (1 - alpha_opt)...
                  .*abs_Y_2(:,lambda);                              %eq.  4
  % compute bias correction
  beta = min(alpha_opt.^2,beta_max);
  mean_P = beta.*mean_P + (1-beta).*P(:,lambda);                    %eq. 20
  mean_P_2 = beta.*mean_P_2 + (1-beta).*P(:,lambda).^2;             %eq. 21
  var_P = mean_P_2 - mean_P.^2;                                     %eq. 22
  inv_Q_eq = var_P./(2*sigma_N_2_prev.^2);                          %eq. 23
  %inv_Q_eq = max(min(inv_Q_eq,INV_Q_EQ_MAX),INV_Q_EQ_MIN/lambda);
  inv_Q_eq = max(min(inv_Q_eq,inv_q_eq_max),1/20/lambda);
  tilda_Qeq_D = (1./inv_Q_eq - 2*M_D)/(1-M_D);                      %eq. 16
  tilda_Qeq_V = (1./inv_Q_eq - 2*M_V)/(1-M_V);                      %eq. 16
  B_min = 1 + 2*(D-1)./tilda_Qeq_D;                                 %eq. 17
  B_min_sub = 1 + 2*(V-1)./tilda_Qeq_V;                             %eq. 17
  % compute inv Q_eq
  inv_Q_eq = mean(inv_Q_eq);
  % minimum search
  B_c = 1 + a_v*sqrt(inv_Q_eq);
  P_x_B_min_x_B_c = P(:,lambda).*B_min*B_c;
  k_mode = P_x_B_min_x_B_c < actmin;
  actmin(k_mode) = P_x_B_min_x_B_c(k_mode);
  actmin_sub(k_mode) = P(k_mode,lambda).*B_min_sub(k_mode)*B_c;
  if subwc == V
    lmin_flag(k_mode) = 0;
    actmin_buff(:,buff_idx) = actmin;
    if buff_idx == U
      buff_idx = 1;
    else
      buff_idx = buff_idx + 1;
    end
    P_min_u = min(actmin_buff,[],2);
    p = find(inv_Q_eq < noise_slope_max_table(:,1));
    if isempty(p)
      noise_slope_max = 1.2;
    else
      noise_slope_max = noise_slope_max_table(p(1),2);
    end
    p = and(and(lmin_flag, actmin_sub < noise_slope_max.*P_min_u),...
            actmin_sub > P_min_u);
    if any(p)
      P_min_u(p) = actmin_sub(p);
      actmin(p) = actmin_sub(p);
    end
    lmin_flag = zeros(L,1);
    subwc = 1;
    actmin(:) = MAX_VALUE;
    actmin_sub(:) = MAX_VALUE;
    sigma_N_2(:,lambda) = P_min_u;
  else
    if subwc > 1
      lmin_flag(k_mode) = 1;
      sigma_N_2(:,lambda) = min(actmin_sub,P_min_u);
      P_min_u = sigma_N_2(:,lambda);
    else
      sigma_N_2(:,lambda) = min(actmin_sub,P_min_u);
    end
    subwc = subwc + 1;
  end
  P_prev = P(:,lambda);
  sigma_N_2_prev = sigma_N_2(:,lambda);
end
end

function [M_D, H_D] = mean_of_the_min_table3(D)
DMH = [  1 0.000 0.000;   2 0.260 0.150;   5 0.480 0.480;   8 0.580 0.780; ...
        10 0.610 0.980;  15 0.668 1.550;  20 0.705 2.000;  30 0.762 2.300; ...
        40 0.800 2.520;  60 0.841 2.900;  80 0.865 3.250; 120 0.890 4.000; ...
       140 0.900 4.100; 160 0.910 4.100];

 p = find(DMH(:,1) == D);
 if ~isempty(p)
   M_D = DMH(p,2);
   H_D = DMH(p,3);
 else 
   p = find(DMH(:,1) < D);
   if isempty(p)
     error('invalid value for: D');
   end
   p = p(end);
   if p > size(DMH,1)-1;
     error('invalid value for: D');
   end
   delta_p_rat = (D - DMH(p,1))/(DMH(p+1,1)-DMH(p,1));
   M_D = DMH(p,2) + (DMH(p+1,2)-DMH(p,2))*delta_p_rat;
   H_D = DMH(p,3) + (DMH(p+1,3)-DMH(p,3))*delta_p_rat;
 end
end
