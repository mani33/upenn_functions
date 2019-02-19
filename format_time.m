function t = format_time(ti)
% t = format_time(neuralynx_time)
% Converts neuralynx time to minutes
ti = double(ti);
t = ((ti-ti(1))*1e-6)/60; % convert microseconds to minutes