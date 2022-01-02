%--------------------------------------------------------------------------
% University of Limerick - Dept. of Electronic and Computer Engineering
%--------------------------------------------------------------------------
% filename:    Coherence_.m
%
% purpose:     This function calculates the frequency Coherence
%              between two signals by applying the equation (3-6)    
%              using the Welch's method as described in section 3.3.2 in Thesis.
%              **It gives the same result to the MATLAB call mscohere().**             
% created by:  Fotios Kostarelos
% created on:  9th March 2021
%--------------------------------------------------------------------------
% Copyright 2021 University of Limerick
%--------------------------------------------------------------------------

function Coherence_ = Coherence_(x,y,z)

% initializations
Fs = 125;
x = x*10^-2;
y = y*10^-2;
x = double(x);
y = double(y);

thetabandStart = 4; 
thetabandEnd = 6; 

betabandStart = 12; 
betabandEnd = 22; 

mybetaCo=zeros(1,betabandEnd-betabandStart);
betaCo=zeros(1,betabandEnd-betabandStart);

mythetaCo=zeros(1,thetabandEnd-thetabandStart);
thetaCo=zeros(1,thetabandEnd-thetabandStart);

Navg = 127;
nfft = 128;
win=hamming(nfft);
window=win.*sqrt(nfft/sum(win.^2));

samplesread = nfft-1;
startingpoint=1;
endpoint=0;
xmulconjyfft = zeros(nfft/2+1,Navg);

magsqmpax= zeros(nfft/2,1);
magsqmpay= zeros(nfft/2,1);

% spectral calculations 
N = length(x);
for i=1:Navg
endpoint=startingpoint+samplesread;
x1=x(startingpoint:endpoint)'.*window;
y1=y(startingpoint:endpoint)'.*window;
xdft1 = fft(x1,nfft);  %<outputs two-sided complex amplitude spectra
ydft1 = fft(y1,nfft);
startingpoint=startingpoint+(samplesread+1)/2;
Sxy1a=ydft1.*conj(xdft1);   
Sxy1c=Sxy1a(1:nfft/2+1);
Sxy1c(2:end-1)=2*Sxy1c(2:end-1);
Sxy1c = Sxy1c';
xmulconjyfft(:,i)=Sxy1c;
end


avgxyfft = mean(xmulconjyfft,2);
avgxyfft= avgxyfft/(nfft*Fs);

startingpointx=1;
endpointx=0;

for i=1:Navg
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

for i= 1:Navg
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
n = 1;
for m=betabandStart:betabandEnd 
 n = n + 1;
 mybetaCo(n)=((real(avgxyfft(m))*real(avgxyfft(m)))+(imag(avgxyfft(m))*imag(avgxyfft(m))))/(P_Hz1x(m)*P_Hz1y(m));
end
summybetaCo = sum(mybetaCo);
mybetaCoh = summybetaCo/(betabandEnd-betabandStart+1);
n = 1;
for m=thetabandStart:thetabandEnd 
 n = n + 1;
 mythetaCo(n)=((real(avgxyfft(m))*real(avgxyfft(m)))+(imag(avgxyfft(m))*imag(avgxyfft(m))))/(P_Hz1x(m)*P_Hz1y(m));
end
summythetaCo = sum(mythetaCo);
mythetaCoh = summythetaCo/(thetabandEnd-thetabandStart+1);

if z==1 
    Coherence_ = mybetaCoh;
else 
    Coherence_ = mythetaCoh;
end 
    

end 
