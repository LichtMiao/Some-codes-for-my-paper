clc;
close all;
clear all;

% network setttings
SbsLambda = 1e-4;
UeLambda = 1e-3;

MbsList = [0,0];
MbsRadius = 200;
SbsRadius = 60;

% system parameters
BW = 150e3;
SbsPtx = 150e-3;
MbsPtx = 750e-3;
Noise_dB = -90;
M = 48; % subcarrier number(M>=UeNum)
B0 = BW/M; % the bandwith of subcarrier
SubcarrierId = (1:M); % 子载波


% channel parameters
alpha = 3;
mean_fading_dB = 0;
std_fading_dB = 49;
PL = @(distance) distance.^(-alpha);

% Cooperation threshold
SinrThreshold_dB = 30;
SinrThreshold = db2lin(SinrThreshold_dB);

% Generate SBS
SbsNum = max([random('poiss',MbsRadius^2.*SbsLambda,1,1), 3]);
SbsLists = random('unif', -MbsRadius, MbsRadius, SbsNum, 2);
% Generate UE
UeNum = random('poiss',MbsRadius^2.*UeLambda,1,1);
UeLists = random('unif', -MbsRadius, MbsRadius, UeNum, 2);

figure;
plot(MbsList(:,1), MbsList(:,2), 'g<'); hold on;
plot(UeLists(:,1), UeLists(:,2), 'r.');  hold on;
plot(SbsLists(:,1), SbsLists(:,2), 'b*');
axis equal; hold off;


% Storage the distance of SBSs and UEs in rectangle
SbsDistance = sqrt(SbsLists(:,1).^2 + SbsLists(:,2).^2);
UeDistance = sqrt(UeLists(:,1).^2 + UeLists(:,2).^2);

% Keep the scope within the macrocell
SbsActiveLists = SbsLists(SbsDistance<MbsRadius, :);
UeActiveLists = UeLists(UeDistance<MbsRadius, :);

figure;
plot(SbsActiveLists(:,1), SbsActiveLists(:,2), 'b*'); hold on;
plot(UeActiveLists(:,1), UeActiveLists(:,2), 'r.');  hold on;
plot(MbsList(:,1), MbsList(:,2), 'g<');
axis equal;

% Caculate the distance from UE to all SBSs
Ue2Sbs = zeros(size(UeActiveLists, 1), size(SbsActiveLists, 1));
for Sbs_iter = 1:size(SbsActiveLists, 1)
    temp = UeActiveLists - repmat(SbsActiveLists(Sbs_iter,:), size(UeActiveLists, 1), 1);
    Ue2Sbs(:,Sbs_iter) = sqrt(temp(:,1).^2 + temp(:,2).^2);
end
% Calculate the distance from UE to the MBS
Ue2Mbs = sqrt(UeActiveLists(:,1).^2 + UeActiveLists(:,2).^2);

Ue_Subcarrier = zeros(size(UeActiveLists,1), M); % 所有用户-所有子载波分配矩阵

[Sue2Sbs, SbsIndex] = min(Ue2Sbs'); % associate with small cell
Sue2Sbs = Sue2Sbs';
SbsIndex = SbsIndex';

UeFlag = zeros(size(UeActiveLists, 1), M);
UeSinr = zeros(size(UeActiveLists, 1), M);
DataRateUe = zeros(size(UeActiveLists, 1), M);

for freq_iter = 1:M
    % Calculate the average Recv power from the MBS
    MbsPrx = (MbsPtx/M).*PL(Ue2Mbs).*FadingCoefficient(db2lin(mean_fading_dB),db2lin(std_fading_dB), size(Ue2Sbs, 1), 1);
    % Calculate the average Recv power from the inter SBS
    SbsInterfPrx = (SbsPtx/M).*PL(Ue2Sbs).*FadingCoefficient(db2lin(mean_fading_dB),db2lin(std_fading_dB), size(Ue2Sbs, 1), size(Ue2Sbs, 2));

    SbsPrx = zeros(size(UeActiveLists, 1), 1);
    InterferenceSbs = zeros(size(UeActiveLists, 1), 1);
    for ue_iter =  1:size(UeActiveLists, 1)
        SbsPrx(ue_iter) = SbsInterfPrx(ue_iter, SbsIndex(ue_iter));
        InterferenceSbs(ue_iter) = sum(SbsInterfPrx(ue_iter,:)) - SbsPrx(ue_iter);
    end
    
    %% 对用户的三种工作模式进行定义,从而对UEs和SBSs进行分类
    UeRatio_dB = lin2db(SbsPrx./MbsPrx); % SINR ratio
   
    for Ue_iter = 1:size(UeActiveLists, 1)
        if UeRatio_dB(Ue_iter)>SinrThreshold_dB
            UeFlag(Ue_iter,freq_iter) = 1; % Sue flag, only associated with the samll cell;
            % If SINR ratio > SinrThreshold, associate to the nearnest SBS
            % Calculate the SINR of SUE
            UeSinr(Ue_iter,freq_iter) = SbsPrx(Ue_iter)./db2lin(Noise_dB);
            DataRateUe(Ue_iter,freq_iter) = B0.*log(1+UeSinr(Ue_iter,freq_iter)).*1e-3;
        elseif (UeRatio_dB(Ue_iter)<SinrThreshold_dB) & (Ue2Sbs(Ue_iter, SbsIndex(Ue_iter))<SbsRadius)
            UeFlag(Ue_iter,freq_iter) = 2; % Cue flag;
            % If SINR ratio < SinrThreshold & Ue2Sbs < SbsRadius, operate on cooperation mode
            % Suffering interference from other SBSs.
            UeSinr(Ue_iter,freq_iter) = (MbsPrx(Ue_iter) + SbsPrx(Ue_iter))/(db2lin(Noise_dB) + InterferenceSbs(Ue_iter,:));
            DataRateUe(Ue_iter,freq_iter) = B0.*log(1+UeSinr(Ue_iter,freq_iter)).*1e-3;
            % Calculate the distance from CUE to MBS
        else
            % Otherwise, associate to MBS
            UeFlag(Ue_iter,freq_iter) = 3; % Mue flag;
            % Calculate the SINR of MUE
            UeSinr(Ue_iter,freq_iter) = MbsPrx(Ue_iter)./(db2lin(Noise_dB));
            DataRateUe(Ue_iter,freq_iter) = B0.*log(1+UeSinr(Ue_iter,freq_iter)).*1e-3;
        end
    end
end

%Subcarrier-UserId allocation
[DateRateFreq, UeId] = max(DataRateUe, [], 1); % 1 maximum value of each column

%% power allocation



