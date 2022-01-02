%--------------------------------------------------------------------------
% University of Limerick - Dept. of Electronic and Computer Engineering
%--------------------------------------------------------------------------
% filename:    Phase_Difference.m
%
% purpose:     This function calculates the frequency Phase difference 
%              between two signals by applying the equation (3-7)    
%              using the Welch's method as described in section 3.3.2 
%              in Thesis.             
% created by:  Fotios Kostarelos
% created on:  9th March 2021
%--------------------------------------------------------------------------
% Copyright 2021 University of Limerick
%--------------------------------------------------------------------------

function Phase_Difference = Phase_Difference(x,y)

% initializations
Fs = 125;

x = x*10^-2;
y = y*10^-2;
N = length(x);
Navg = 127;
nfft = 128;
win=hamming(nfft);
window=win.*sqrt(nfft/sum(win.^2));

betabandStart = 12; 
betabandEnd = 22;

betaPD1=zeros(1,betabandEnd-betabandStart);

samplesread = nfft-1;
startingpoint=1;
endpoint=0;
xmulconjyfft = zeros(nfft/2+1,Navg);

% spectral calculations
for i=1:Navg
endpoint=startingpoint+samplesread;
x1=x(startingpoint:endpoint)'.*window;
y1=y(startingpoint:endpoint)'.*window;
xdft1 = fft(x1,nfft);  %<outputs two-sided complex amplitude spectra
ydft1 = fft(y1,nfft);
startingpoint=startingpoint+(samplesread+1)/2;
Sxy1a=ydft1.*conj(xdft1);   
Sxy1b=Sxy1a./N^2;
Sxy1c=Sxy1b(1:nfft/2+1);
Sxy1c(2:end-1)=2*Sxy1c(2:end-1);
Sxy1c = Sxy1c';
xmulconjyfft(:,i)=Sxy1c;
end

avgxyfft = mean(xmulconjyfft,2);
avgxyfft= avgxyfft*((length(x).^2)/(nfft*Fs));


% feature calculations
n=1;
for m=betabandStart:betabandEnd 
 n = n + 1;
 betaPD1(n)=atan(imag(avgxyfft(m))/real(avgxyfft(m)))/19; %% Reference thatcher's paper & theory
end
sumbetaPD1 = sum(betaPD1);

Phase_Difference = sumbetaPD1/(betabandEnd-betabandStart+1);

end
