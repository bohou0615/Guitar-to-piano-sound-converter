%% EE6641 HW: Linear prediction and Levinson-Durbin k-parameter estimation
% Created May 2013.
% Last updated May 8, 2019 for this year's HW5.
% Yi-Wen Liu
clear; close all;

DIR = './';
DIR_guitar = 'music_guitar/';
DIR_piano = 'music_piano/';
%FILENAME = 'u2.wav';
%FILENAME = 'a-male-singing.wav';
%FILENAME = 'i-male-singing.wav';
%FILENAME = 'a_1.wav';
%FILENAME1 = 'piano.wav';
%FILENAME2 = 'guitar.wav';
FILENAME1 = 'guitar_c4.wav';
FILENAME2 = 'piano_c4.wav';

%% Parameters to play with
framelen = 1.5; % second. [INVESTIGATE]
p = 200; % linear prediction order. [INVESTIGATE]

%% 1
%[y1,fs1] = audioread([DIR FILENAME1]);
[y1,fs1] = audioread([DIR FILENAME1]);
if size(y1,2) == 2
    y1 = (y1(:,1) + y1(:,2))/2;
end
%soundsc(y,fs1);
fs = 16000;
y1 = resample(y1,fs,fs1);
L = framelen*fs;

if L<=p
    disp('Linear prediction requires the num of equations to be greater than the number of variables.');
end

sw.emphasis = 1; % default = 1

numFrames1 = floor(length(y1)/L);

%% 2
%[y2,fs2] = audioread([DIR FILENAME2]);
[y2,fs2] = audioread([DIR FILENAME2]);
if size(y2,2) == 2
    y2 = (y2(:,1) + y2(:,2))/2;
end

y2 = resample(y2,fs,fs2);
L = framelen*fs;

if L<=p
    disp('Linear prediction requires the num of equations to be greater than the number of variables.');
end

sw.emphasis = 1; % default = 1
numFrames2 = floor(length(y2)/L);

numFrames = min(numFrames1,numFrames2);

%% start
excitat1 = zeros(size(y1));
e_n1 = zeros(p+L,1);
yyy1 = zeros(size(y1));
y_n1 = zeros(p+L,1);


LPcoeffs1 = zeros(p+1,numFrames);
Kcoeffs1 = zeros(p,numFrames); % reflection coeffs

Nfreqs = 2^nextpow2(2*L-1)/2; % Num points for plotting the inverse filter response
df = fs/2/Nfreqs;
ff = 0:df:fs/2-df;

if sw.emphasis == 1,
    y_emph1 = filter([1 -0.95],1,y1); 
                %[PARAM] -0.95 may be tuned anywhere from 0.9 to 0.99
else
    y_emph1 = y1;
end

%% Linear prediction and estimation of the source e_n
win1 = ones(L,1); H_peak1 = zeros(1,numFrames);% Rectangular window.


for kk = 1:numFrames
    ind = (kk-1)*L+1:kk*L;
    ywin1 = y_emph1(ind).*win1;
    Y1 = fft(ywin1,2*Nfreqs);
    R1 = ifft(Y1.*conj(Y1));
 
    %% [INVESTIGATE] You are encouraged to try implementing the
    % Levinson-Durbin Algorithm, and check your results against levinson().
    % ...
   % A1 = lpc(ywin1,p);
    % Coding your "Mylevinson()", you can compute the LP coefficients, error and
    % reflection coefficients.
    [A1,errvar1,K1] = levinson(R1,p); %(Matlab built-in function for levinson-durbin algorithm.)
    %[A,errvar,K] = Mylevinson(R,p);     
    %% Preparation for data visualization
    if kk == 1,      
        e_n1(p+1:end) = filter(A1,[1],ywin1);
    else
        ywin_extended = y1((kk-1)*L+1-p:kk*L);  %% WORTH TWEAKING
        e_n1 = filter(A1,[1],ywin_extended);
    end
    excitat1(ind) = e_n1(p+1:end);
%   soundsc(excitat1,fs); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ewin1 = excitat1(ind).*win1;
    if kk == 1      
        y_n1(p+1:end) = filter(1,A1,ewin1);
    else
        y_n1 = filter(1,A1,e_n1);
    end
    yyy1(ind) = y_n1(p+1:end);
    
    LPcoeffs1(:,kk) = A1;
    %Kcoeffs1(:,kk) = K1;
end
%soundsc(yyy1,fs);
%audiowrite('output.wav',yyy1, fs);
audiowrite('excitat_guitar.wav',excitat1, fs);



%% start 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excitat2 = zeros(size(y2));
e_n2 = zeros(p+L,1);
yyy2 = zeros(size(y2));
y_n2 = zeros(p+L,1);

LPcoeffs2 = zeros(p+1,numFrames);
Kcoeffs2 = zeros(p,numFrames); % reflection coeffs

if sw.emphasis == 1,
    y_emph2 = filter([1 -0.95],1,y2); 
                %[PARAM] -0.95 may be tuned anywhere from 0.9 to 0.99
else
    y_emph2 = y2;
end

%% Linear prediction and estimation of the source e_n
win2 = ones(L,1); H_peak2 = zeros(1,numFrames);% Rectangular window.
%{
LPcoeffs1(:,1) = LPcoeffs1(:,4);
LPcoeffs1(:,2) = LPcoeffs1(:,4);
LPcoeffs1(:,3) = LPcoeffs1(:,4);
LPcoeffs1(:,5) = LPcoeffs1(:,4);
LPcoeffs1(:,6) = LPcoeffs1(:,4);
%}
for kk = 1:numFrames
    ind = (kk-1)*L+1:kk*L;
    ywin2 = y_emph2(ind).*win2;

    Y2 = fft(ywin2,2*Nfreqs);
    R2 = ifft(Y2.*conj(Y2));
    
   % A2 = lpc(ywin2,p);
  
    [A2,errvar2,K2] = levinson(R2,p); %(Matlab built-in function for levinson-durbin algorithm.)

%% Preparation for data visualization
    if kk == 1,      
        e_n2(p+1:end) = filter(A2,[1],ywin2);
    else
        ywin_extended = y2((kk-1)*L+1-p:kk*L);  %% WORTH TWEAKING
        e_n2 = filter(A2,[1],ywin_extended);
    end
    excitat2(ind) = e_n2(p+1:end);
    
    
    ewin2 = excitat2(ind).*win2;
    if kk == 1
        y_n2(p+1:end) = filter(1,LPcoeffs1(:,kk),ewin2);%%%%%%%%%%%%%%%%%%%
    else
        y_n2 = filter(1,LPcoeffs1(:,kk),e_n2);          %%%%%%%%%%%%%%%%%%%
    end
    yyy2(ind) = y_n2(p+1:end);
    
    LPcoeffs2(:,kk) = A2;
   % Kcoeffs2(:,kk) = K2;
end
%soundsc(excitat1,fs);
%soundsc(excitat2,fs);
%{
figure(1)
subplot(211);
plot((excitat1));
subplot(212);
plot((excitat2));
figure(2)
subplot(211);
plot((y1));
subplot(212);
plot((y2));
%}
audiowrite('excitat_piano.wav',excitat2, fs);
audiowrite('output_gtop.wav',yyy2,fs);%guitar to piano
%% [INVESTIGATE] play the estimated source signal. 
% With a proper choice of LP order, frame length, and the pre-emphasis filter coefficient, 
% the sound should lack the original vowel quality in x[n].
%soundsc(excitat2,fs); 

soundsc(yyy2,fs); 
