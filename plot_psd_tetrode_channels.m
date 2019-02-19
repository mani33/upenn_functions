function plot_psd_tetrode_channels(dataFolder,xlimit)

if nargin<2
    xlimit = [0 16];
end
f = dir(fullfile(dataFolder,'t*c*.ncs'));
by = [f.bytes];
s = by>16384;
f = f(s);
fn = {f.name};
nf = length(fn);
for i = 1:nf
    br = baseReaderNeuralynx(fullfile(dataFolder,fn{i}));
    x = br(1:100000,1);
    pwelch(x,[],[],[],32000);
    xlim(xlimit)
    pause
end;
    