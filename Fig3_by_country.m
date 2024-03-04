% This function replicates Figure 3 from the article
%
% @article{TAVARESGARCIA2024105006,
% title = {The impact of monetary policy on income inequality: Does inflation targeting matter?},
% author = {Francisco {Tavares Garcia} and Jamie L. Cross},
% journal = {Finance Research Letters},
% doi = {https://doi.org/10.1016/j.frl.2024.105006},
% }
%
% And it has minor modifications to functions from the article
%
% @article{MUMTAZ2017410,
% title = {The impact of monetary policy on inequality in the UK. An empirical analysis},
% author = {Haroon Mumtaz and Angeliki Theophilopoulou},
% journal = {European Economic Review},
% doi = {https://doi.org/10.1016/j.euroecorev.2017.07.008},
% }

clear;clc;
addpath(genpath('.\'));
sheets = sheetnames('\data\fullset.xlsx');

%%% SELECT COUNTRY %%%
% Generating all IRFs from Figure 3 can take a couple of days on a personal
% computer. So, we split the analysis by country, which can be selected
% from the list below:

% country = "ca"; % Canada
% country = "fr"; % France
% country = "ge"; % Germany
% country = "it"; % Italy
% country = "jp"; % Japan
country = "uk"; % United Kingdom
% country = "us"; % United States

tablenames = [];
for i = 1:size(sheets, 1)
    ct = strsplit(sheets(i),' ');
    if ct(1) == country
        tablenames = [tablenames; sheets(i)];
    end
end

% Loop to analyse each sheet of the income.xlsx spreadsheet.
% It allows for multiple dataset IRFs in a single run.
for i = 1:size(tablenames, 1)
    table = readtable('\data\fullset.xlsx', 'Sheet', tablenames(i));

    % inputs
    % names = table.Properties.VariableNames(2:6);    % column names
    names = {'GDP PER CAPITA', 'CPI', 'GINI', 'TBILL', 'NEER'};
    data = table.Variables; % data variable
    data(:,1) = [];         % remove variable dates

    reps=15000; %total reps, length: over 12 minutes (6 Threads)
    % reps=1500;
    burn=10000; %burn in
    % burn=1000;
    % update=1;%prints every update iter
    update=reps/100; %prints every 100 update iter

    maxtrys=1000; %max tries for stable draw and to find A0 matrix
    horizon=40; %forecast horizon
    L=4;        % number of lags in the VAR

    % For the USA and Japan, 2 lags were used due to the small sample post
    % inflation targeting adoption
    if (country == "jp") || (country == "us")
        L=2;
    end

    identification=2; %1 for Cholesky ordering as in data, any other number for sign restrictions
    %priors Banbura et al. JAE 2009
    lamdaP=0;% tightness of prior on lags
    tauP=10*lamdaP; %tightness of prior on sum of coefficients
    epsilonP=1/1000; %tightness of prior on constant
    mreps=1000;

    %If identification is ~=1 specify sign restrictions
    pattern=zeros(5,5);

    signs=[-1 -1 0  1 1 ];
    pattern(4,:)=signs;
    timemat=zeros(rows(pattern),cols(pattern));


    disp(tablenames(i))
    %this function estimates a VAR with gibbs sampling. For draw after burn it
    %tries to find 1 A0 matrix that satisfies sign restriction. Other
    %approaches are described in
    %http://www-personal.umich.edu/~lkilian/km080210.pdf
    tic
    [fsave,hsave,bsave,fvsave,emat]=bvar_sign(data,pattern,timemat,L,reps,burn,horizon,update,maxtrys,lamdaP,tauP,epsilonP,identification,mreps,tablenames(i));
    toc

    % plot the response to shock 5 the policy shock
    f = figure('Name', tablenames(i));
    temp=squeeze(fsave(:,4,:,:));
    tt=0:horizon-1;
    for j=1:cols(data)
        subplot(3,2,j)
        temp1=squeeze(prctile(temp(:,:,j),[50  16 84],1)); % 68% CI
        % temp1=squeeze(prctile(temp(:,:,j),[50  2.5 97.5],1)); % 95% CI
        plotx2(tt,temp1')
        title(names{j});
    end

    % saves figure
    figurename = strcat('.\results\irf\', tablenames(i), '_', string(L), 'lags.fig');
    % figurename = strcat('.\results\robustness\', tablenames(i), '_', string(L), 'lags.fig');
    saveas(f, figurename);
end
