function pstr = get_plessthan_str(p)
% pstr = get_plessthan_str(p)
% gives p < 0.01 like ones as output

if p < 0.00001
    pstr = 'p < 0.00001';
elseif p < 0.0001
    pstr = 'p < 0.0001';
elseif p < 0.001
    pstr = 'p < 0.001';
elseif p < 0.01
    pstr = 'p < 0.01';
else
    pstr = sprintf('p = %0.3f',p);
end
disp('----------------------------------')
fprintf('Actual p-value is: %0.10f\n',p)
disp('----------------------------------')