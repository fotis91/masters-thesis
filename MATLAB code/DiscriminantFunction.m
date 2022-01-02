%--------------------------------------------------------------------------
% University of Limerick - Dept. of Electronic and Computer Engineering
%--------------------------------------------------------------------------
% filename:    Discriminant_Function.m
%
% purpose:     Discriminant Function realisation as described in thesis. 
%              Uses the MATLAB functions which are the features for 
%              Relative_Power, Amplitude_Asymmetry, Coherence_ and 
%              Phase_Difference, to calculate the 20 variables from 
%              Fig 4-13.
%              
% created by:  Fotios Kostarelos
% created on:  9th March 2021
%--------------------------------------------------------------------------
% Copyright 2021 University of Limerick
%--------------------------------------------------------------------------

%clear all
%close all
%clc

N=2^13; % input length 

% theta coherence channels Fp1 F3
x = EEG.data(2,N+1:2*N);
y = EEG.data(2,1:N);
Theta_Coherence_between_Fp1_F3 = Coherence_(x,y,2);

% beta coherence channels T7 P7
x = EEG.data(3,N+1:2*N);
y = EEG.data(1,1:N);
Beta_Coherence_between_T7_P7 = Coherence_(x,y,1);

% beta coherence channels C3 P3
x = EEG.data(4,N+1:2*N);
y = EEG.data(4,1:N);
Beta_Coherence_between_C3_P3 = Coherence_(x,y,1);

% beta phase difference channels Fp2 F4
x = EEG.data(5,N+1:2*N);
y = EEG.data(6,N+1:2*N);
Beta_Phase_Difference_between_Fp2_F4 = Phase_Difference(x,y);

% beta phase difference channels F3 F4
x = EEG.data(2,1:N);
y = EEG.data(6,N+1:2*N);
Beta_Phase_Difference_between_F3_F4 = Phase_Difference(x,y);

% alpha,beta amplitude asymmetry channels F4 P8
x = EEG.data(6,N+1:2*N);
y = EEG.data(3,1:N);
Alpha_Amplitude_Asymmetry_between_F4_P8 = Amplitude_Asymmetry(x,y,1);
Beta_Amplitude_Asymmetry_between_F4_P8 = Amplitude_Asymmetry(x,y,2);

% alpha,beta amplitude asymmetry channels F8 P8
x = EEG.data(8,N+1:2*N);
y = EEG.data(3,1:N);
Alpha_Amplitude_Asymmetry_between_F8_P8 = Amplitude_Asymmetry(x,y,1);
Beta_Amplitude_Asymmetry_between_F8_P8 = Amplitude_Asymmetry(x,y,2);

% alpha,beta amplitude asymmetry channels F4 O2
x = EEG.data(6,N+1:2*N);
y = EEG.data(6,1:N);
Alpha_Amplitude_Asymmetry_between_F4_O2 = Amplitude_Asymmetry(x,y,1);
Beta_Amplitude_Asymmetry_between_F4_O2 = Amplitude_Asymmetry(x,y,2);

% alpha amplitude asymmetry channels F3 O1
x = EEG.data(2,1:N);
y = EEG.data(5,1:N);
Alpha_Amplitude_Asymmetry_between_F3_O1 = Amplitude_Asymmetry(x,y,1);

% alpha amplitude asymmetry channels O1 F7
x = EEG.data(5,1:N);
y = EEG.data(1,N+1:2*N);
Alpha_Amplitude_Asymmetry_between_O1_F7 = Amplitude_Asymmetry(x,y,1);

% alpha relative power channel P3
x = EEG.data(4,1:N);
Alpha_Relative_Power_for_P3 = Relative_Power(x);

% alpha relative power channel P4
x = EEG.data(7,1:N);
Alpha_Relative_Power_for_P4 = Relative_Power(x);

% alpha relative power channel O1
x = EEG.data(5,1:N);
Alpha_Relative_Power_for_O1 = Relative_Power(x);

% alpha relative power channel O2
x = EEG.data(6,1:N);
Alpha_Relative_Power_for_O2 = Relative_Power(x);

% alpha relative power channel T8
x = EEG.data(8,1:N);
Alpha_Relative_Power_for_T8 = Relative_Power(x);

% alpha relative power channel P7
x = EEG.data(1,1:N);
Alpha_Relative_Power_for_P7 = Relative_Power(x);

% alpha relative power channel P8
x = EEG.data(3,1:N);
Alpha_Relative_Power_for_P8 = Relative_Power(x);
