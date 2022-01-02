%--------------------------------------------------------------------------
% University of Limerick - Dept. of Electronic and Computer Engineering
%--------------------------------------------------------------------------
% filename:    Amplitude_Asymmetry.m
%
% purpose:     This function calculates the frequency Amplitude Asymmetry
%              between two signals by applying the equation (3-5)    
%              using the Welch's method as described in section 3.3.2 
%              in Thesis.             
% created by:  Fotios Kostarelos
% created on:  9th March 2021
%--------------------------------------------------------------------------
% Copyright 2021 University of Limerick
%--------------------------------------------------------------------------

function Amplitude_Asymmetry = Amplitude_Asymmetry(x,y,z)

% variables initializations 
Fs = 125; 
Navg = 127;
nfft = 128;
samplesread = nfft-1;
alphabandStart = 7; 
alphabandEnd = 13; 
betabandStart = 12; 
betabandEnd = 22; 
startingpointx=1;
endpointx=0;

% vectors initializations
win=hamming(nfft);
window=win.*sqrt(nfft/sum(win.^2));

magsqmpax= zeros(nfft/2,1);
magsqmpay= zeros(nfft/2,1);

betaAmplitudeAsymmetry=zeros(1,betabandEnd-betabandStart);
alphaAmplitudeAsymmetry=zeros(1,alphabandEnd-alphabandStart);

% scaling(not part of the HLS code)
x = double(x);
y = double(y);
x = x*10^-2;
y = y*10^-2;


% spectral calculations 
for k=1:Navg
    endpointx=startingpointx+samplesread;
    x_column1=x(startingpointx:endpointx)';
    startingpointx=startingpointx+(samplesread+1)/2;
    xw = x_column1.*window;               % apply window of length nfft
    X= fft(xw);                           % DFT
    X= X(1:nfft/2);                      % retain samples from 0 to fs/2
    magsqx= real(X).^2 + imag(X).^2;      % DFT magnitude squared
    magsqmpax(:,1) = magsqmpax(:,1) + magsqx;% sum of DFT mag squared
end

magsqmpax = magsqmpax/Navg;                 % average of DFT mag squared
P_bin1x= 2/nfft.^2 *magsqmpax;              % W/bin power spectrum                     
P_Hz1x= P_bin1x*nfft/Fs;                    % W/Hz power spectrum
P_Hz1x(1)=P_Hz1x(1)/2;    

startingpointy=1;
endpointy=0;

for l= 1:Navg
    endpointy=startingpointy+samplesread;
    y_column1=y(startingpointy:endpointy)';
    startingpointy=startingpointy+(samplesread+1)/2;
    yw = y_column1.*window;               % apply window of length nfft
    Y= fft(yw);                           % DFT
    Y= Y(1:nfft/2);                      % retain samples from 0 to fs/2
    magsqy= real(Y).^2 + imag(Y).^2;      % DFT magnitude squared
    magsqmpay(:,1) = magsqmpay(:,1) + magsqy;% sum of DFT mag squared
end
            
magsqmpay = magsqmpay/Navg;                 % average of DFT mag squared
P_bin1y= 2/nfft.^2 *magsqmpay;              % W/bin power spectrum                     
P_Hz1y= P_bin1y*nfft/Fs;                    % W/Hz power spectrum
P_Hz1y(1)=P_Hz1y(1)/2;

% feature calculations
n=1;
for m=betabandStart:betabandEnd 
 n = n + 1;
 betaAmplitudeAsymmetry(n)=(sqrt(P_Hz1x(m))-sqrt(P_Hz1y(m)))/(sqrt(P_Hz1x(m))+sqrt(P_Hz1y(m))); 
end
sumbetaAA = sum(betaAmplitudeAsymmetry);
betaAAWelch = sumbetaAA/(betabandEnd-betabandStart+1);


n = 1;
for m=alphabandStart:alphabandEnd 
 n = n + 1;
 alphaAmplitudeAsymmetry(n)=(sqrt(P_Hz1x(m))-sqrt(P_Hz1y(m)))/(sqrt(P_Hz1x(m))+sqrt(P_Hz1y(m))); 
end
sumalphaAA = sum(alphaAmplitudeAsymmetry);
alphaaAAWelch = sumalphaAA/(alphabandEnd-alphabandStart+1);

% output 
if z==1
    Amplitude_Asymmetry = alphaaAAWelch;
else 
    Amplitude_Asymmetry = betaAAWelch;
end 

end 
