function varargout = plotfft(x,Fs)

L = length(x);
hL = round(L/2);
Mr = 2/L; % Magnitude adjustment factor (see Lyon's DSP book)
w = Mr * abs(fft(x));
fresol = Fs/L;
f = (0:hL-1)*fresol;
vw = (w(1:hL));
plot(f,vw,'.-')
xlabel('Freq')
ylabel('Amplitude')
if nargout==1
    varargout{1} = vw;
end
if nargout==2
    varargout{1} = vw;
    varargout{2} = f;
end
    