function [wps,Fc_vec] = get_wavelet_ps(data_vec,Fs,varargin)
% [wps,t,Fc_vec] = get_wavelet_ps(data_vec,Fs,df)
% Mani Subramaniyan
% 2019-03-08
args.freq_std_scaling_meth = 'log'; % can be 'constant','log','linear'
args.freq_std_bound = [0.5 5]; % beginning and end of std of frequency response
args.freq_scaling_meth = 'log'; % can be 'log', 'linear'
args.freq_bound = [1 80];
args.nFreq = 80;
args.wave_width_t = 1.5;
args.plot_wavelet = 0;
args.meth = 'conv';
args.zscore = 1;
args.zscoredim = 1; % zscore along rows or columns
args = parseVarArgs(args,varargin{:});

data_vec = double(data_vec);
N = length(data_vec);
% Determine the bandwidth param
Fc_vec = get_Fc(args);
Fb_vec = get_Fb(args);
% Fc_vec = 10*ones(1,length(Fc_vec));
wave_nsamples = 2*round(args.wave_width_t*Fs)+1; % make sure it is odd number
wps = zeros(args.nFreq,N);

%% Prepare for FFT based method
nConv = N + wave_nsamples -1;
halfWavelet = (wave_nsamples-1)/2;
fft_data = fft(data_vec,nConv);

for i = 1:args.nFreq
    % Create wavelet
    psi = get_wavelet(-args.wave_width_t,args.wave_width_t,wave_nsamples,Fb_vec(i),Fc_vec(i));
    % Make sure mean that norm is 1
    psi = psi/norm(psi);
    fft_psi = fft(psi,nConv);    
    % Get time domain signal back after ifft
    conv_out = abs(ifft(fft_data.*fft_psi,nConv));  
    conv_out = conv_out(halfWavelet+1:end-halfWavelet); 
    if args.plot_wavelet
        F = (0:nConv-1)*Fs/nConv;
        plot(F,abs(fft_psi));
        xlim([0 args.freq_bound(end)])
        pause
    end
    switch args.meth
        case 'conv'
            % Convolve data
            wps(i,:) = abs(conv(data_vec,psi,'same'));
        case 'fft'
            wps(i,:) = conv_out;
    end
end
if args.zscore
    wps = zscore(wps,[],args.zscoredim);
end


function Fc_vec = get_Fc(args)
fbound = args.freq_bound;
switch args.freq_scaling_meth
    case 'log'
        Fc_vec = logspace(log10(fbound(1)),log10(fbound(2)),args.nFreq);
    case 'linear'
        Fc_vec = linspace(fbound(1),fbound(2),args.nFreq);
    otherwise
        error('Unknown method')
end
function Fb_vec = get_Fb(args)
fbound = args.freq_std_bound;
switch args.freq_std_scaling_meth
    case 'constant'
        Fb_min = 1/(2*pi*fbound(1)^2);
        Fb_max = Fb_min;
        Fb_vec = linspace(Fb_min,Fb_max,args.nFreq);
    case 'log'
        st = logspace(log10(fbound(1)),log10(fbound(2)),args.nFreq);
        Fb_vec = 1./(2*pi*st.^2);
    case 'linear'
        st = linspace(fbound(1),fbound(2),args.nFreq);
        Fb_vec = 1./(2*pi*st.^2);
    otherwise
        error('Unknown method')
end

%%
function [wt,t] = get_wavelet(LB,UB,N,Fb,Fc)

t = linspace(LB,UB,N)';
wt = (1/sqrt(pi*Fb))*exp((1i*2*pi*Fc*t)-((t.^2)/Fb));

