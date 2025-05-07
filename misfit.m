function [misfit] = misfit(rho_obs, phase_obs, rho_cal, phase_cal)
%MISFIT Calculates the average misfit between observed and calculated MT data
%
% INPUT:
%   rho_obs   : Vector of observed apparent resistivity values
%   phase_obs : Vector of observed phase values (in degrees)
%   rho_cal   : Vector of calculated apparent resistivity values (from model)
%   phase_cal : Vector of calculated phase values (from model)
%
% OUTPUT:
%   misfit : Average combined misfit

% Convert degrees to radians (for phase calculation)
d2r = pi / 180;

% Number of data points
nd = length(rho_obs);

% Initialize misfit vector
m = zeros(1, nd);

% Loop through each data point
for i = 1:nd
    % Compute misfit as the sum of:
    % 1. Absolute logarithmic difference in resistivity (log10 scale)
    % 2. Absolute difference in phase (in radians)
    m(i) = abs(log10(rho_cal(i) / rho_obs(i))) + ...
           abs(d2r * phase_cal(i) - d2r * phase_obs(i));
end

% Final misfit is the average over all data points
misfit = sum(m) / nd;

end
