%--------------------------------------------------------------------------
% University of Limerick - Dept. of Electronic and Computer Engineering
%--------------------------------------------------------------------------
% filename:    Relative_Power.m
%
% purpose:     This function calculates the frequency Relative Power 
%              for a signal by applying the equation (3-4)    
%              using the Welch's method as described in section 3.3.2.
%              in Thesis.
%                           
% created by:  Fotios Kostarelos
% created on:  9th March 2021
%--------------------------------------------------------------------------
% Copyright 2021 University of Limerick
%--------------------------------------------------------------------------

function Relative_Power = Relative_Power(x)

% initializations 
Fs=125; %% sampling rate

Navg = 127; %% Welch's method iterations for averaging
nfft = 128; %% FFT points
win=hamming(nfft); %% window function
window=win.*sqrt(nfft/sum(win.^2));
samplesread = nfft-1;

magsqmpax= zeros(nfft/2,1);


alphabandStart = 7; 
alphabandEnd = 13; 

TotalbandStart = 1;
TotalbandEnd = 22;


x = x*10^-2;
startingpointx=1;
endpointx=0;

% spectral calculations 
for i=1:Navg 
    endpointx=startingpointx+samplesread;
    x_column1=x(startingpointx:endpointx)'; %% Reading data from input equal in size with FFT points
    startingpointx=startingpointx+(samplesread+1)/2;
    xw = x_column1.*window;               % apply window of length nfft
    X= fft(xw,nfft);                           % DFT
    X= X(1:nfft/2);                      % retain samples from 0 to fs/2
    magsqx= real(X).^2 + imag(X).^2;      % DFT magnitude squared
    magsqmpax(:,1) = magsqmpax(:,1) + magsqx;% sum of DFT mag squared
end

magsqmpax = magsqmpax/Navg;                 % average of DFT mag squared
P_bin1x= 2/nfft.^2 *magsqmpax;              % W/bin power spectrum                     
P_Hz1x= P_bin1x*nfft/Fs;                    % W/Hz power spectrum
P_Hz1x(1)=P_Hz1x(1)/2;         

% feature calculations
overallPower = sum(P_Hz1x(TotalbandStart:TotalbandEnd));
alphaPower = sum(P_Hz1x(alphabandStart:alphabandEnd));

Relative_Power = (alphaPower/overallPower) * 100;

end
