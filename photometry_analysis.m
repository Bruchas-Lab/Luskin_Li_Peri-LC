% This code will give you event-related photometry graphs and the basic
% stats for pre/post event (at the end)

%%

clear all
close all
clc

%%

experimentlabel = 'VGAT food deprived feeding stop';

k = 1;

for sess = [1:NUMBER OF SESSIONS] % fill in number of sessions

    clearvars -except sess k itchStack experimentlabel
    
% Put session IDs here

if sess == 1
photoname = 'MOUSENAME_1.mat'; % extracted photometry file
pokeTimes1 = dlmread('MOUSENAME_1_EATING_TIMES.txt'); % behavior stored in .txt file with timestamps in first column and "1" after a tab
end
if sess == 2
photoname = 'MOUSENAME_2.mat';
pokeTimes1 = dlmread('MOUSENAME_2_EATING_TIMES.txt');
end
% list out rest of sessions here

%%


load(photoname)

%%  

pokeTimes = round(pokeTimes1(:,1));
eventtype = round(pokeTimes1(:,2));
pokeTimes = pokeTimes(eventtype==1);  

experDuration = floor(max(Dts));

%%

 

% get photometry trace
Fs = 1017.25;

C1 = zeros(experDuration,1);
C1(pokeTimes) = 1;  % indicate Score times

stopdelay = 0;
load(photoname)
photom1 = data1(round(Fs*stopdelay)+1:round(Fs*(experDuration+stopdelay)));
time = linspace(1/Fs,experDuration,experDuration*Fs);

behavior = C1;
timeBehav = linspace(0,experDuration,length(behavior));
 
    % Z score photometry data
   photom1 = (photom1-mean(photom1))./std(photom1);


    
    figure(sess+100)
    
    plot(timeBehav,behavior+10,'g',...
        decimate(time,1000),decimate(photom1,1000));
    
    ylabel('Z score                       Events per second')
    xlabel('Time (sec)')
    yticks([-10 -5 0 5 10 15 20 25 30 35 40])
    yticklabels({'-10','-5','0','5','0','5','10','15','20','25','30'})
    
    windtop = 20;
     axis([0 experDuration -5 windtop])



%%

n = 1;
for p = 2:length(behavior)

    if behavior(p-1) == 0
        if behavior(p) == 1
        
           licks1(n) = timeBehav(p);
           n = n + 1; 
           
        end
    end
    
end

  
  timewindow = 30; % +/- seconds from event

%% triggered averaging (of photom trace)

decifactor = 1;

photom1 = decimate(photom1,decifactor);
time = decimate(time,decifactor);
% 
Fs = 1017.25/decifactor;


samplewindow = round(timewindow*Fs); % samples per time period

licks2 = licks1(licks1<(experDuration-timewindow));
licks = licks2(licks2>timewindow);


LickTrig = zeros(length(licks),2*samplewindow);

for p = 1:length(licks)
    
    trigindex = zeros(1,length(photom1));
    
    [~,stopidx] = min(abs(time-(licks(p)-timewindow)));
    trigindex(stopidx:(stopidx+size(LickTrig,2)-1)) = 1;
    
    LickTrig(p,:) = photom1(trigindex==1);
    
end

photoPerLick = mean(LickTrig,1);
%

sem = std(LickTrig,0,1)./sqrt(size(LickTrig,1)); % sem = std/sqrt(n)

% event triggered average
figure(sess+200)

hold on
x = decimate(linspace(-timewindow,timewindow,length(photoPerLick)),300);
y = decimate(photoPerLick,300);
eb = decimate(sem,300);
lineProps.col{1} = 'blue';
mseb(x,y,eb,lineProps,1);
L = line([0 0],[-3 3]);
set(L,'Color','black')
xlabel('Peri-Event Time (sec)')
ylabel('Z score (smoothed)')
hold off

%% heat map (x dim: time, ydim: event, zdim: deltaF/F)
figure(sess+300)

hold on
imagesc(linspace(-timewindow,timewindow,length(photoPerLick)),1:size(LickTrig,1),LickTrig)
L = line([0 0],[0 length(licks)+1]);
set(L,'Color','black')
xlabel('Peri-Event Time (sec)')
ylabel('Bout Number')
cb = colorbar;
title(cb,'Z score')
caxis([-1.5 1.5])
colormap(flipud(brewermap([],'YlGnBu')))
xlim([-timewindow timewindow])
ylim([0 length(licks)+1])
hold off

%%

photom10hz = resample(photom1,1,floor(length(photom1)/length(behavior)));

photom10hz = photom10hz(1:length(timeBehav));

baseAVG = mean(photom10hz(behavior~=1)); 
eatAVG = mean(photom10hz(behavior==1));

baseSEM = std(photom10hz(behavior~=1))/sqrt(length(photom10hz(behavior~=1)));
eatSEM = std(photom10hz(behavior==1))/sqrt(length(photom10hz(behavior==1)));

[~,p] = ttest2(photom10hz(behavior~=1),photom10hz(behavior==1));

itchStack{k} = LickTrig;

k = k + 1;

end

% close all

%% combined plot

superStack = [];

    for  n = 1:length(itchStack)
    superStack = cat(1,superStack,itchStack{n});
    end

LickTrig1 = superStack;



photoPerLick = mean(LickTrig1,1);
%

sem = std(LickTrig1,0,1)./sqrt(size(LickTrig1,1)); % sem = std/sqrt(n)

% event triggered average
close all
figure(61)

hold on
x = decimate(linspace(-timewindow,timewindow,length(photoPerLick)),1000);
y = decimate(photoPerLick,1000);
eb = decimate(sem,1000);
lineProps.col{1} = [0 0.5 0];
mseb(x,y,eb,lineProps,1);
L = line([0 0],[-3 5]);
axis([-10 20 -0.6 0.9])
yticks([-0.4 0 0.4 0.8])
xticks([-10 0 10 20])
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('Activity (z scored)')
title(experimentlabel)
hold off

%print -painters -depsc file_name_here.eps %quick dirty export


%% heat map (x dim: time, ydim: event, zdim: deltaF/F)
figure(62)

hold on
imagesc(linspace(-timewindow,timewindow,length(photoPerLick)),1:size(LickTrig1,1),LickTrig1)
L = line([0 0],[0 size(LickTrig1,1)+1]);
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('Bout Number')
title(experimentlabel)
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
xlim([-20 60])
ylim([0 size(LickTrig1,1)+1])
hold off

%print -painters -depsc heatmapSuM.eps

%%
%% pre versus post event statistics

%%%%%%%%%% Pick time window
preWindow = [-10 0];
postWindow = [0 10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = 1017.25;

preIndex = round(preWindow.*Fs+size(LickTrig1,2)/2);
postIndex = round(postWindow.*Fs+size(LickTrig1,2)/2);

% slice LickTrig1 2d matrix and take mean of each interval
preList = mean(LickTrig1(:,preIndex(1):preIndex(2)),2);
postList = mean(LickTrig1(:,postIndex(1):postIndex(2)),2);

preList
postList

preMean = mean(preList)
preSTD = std(preList)
preN = length(preList)

postMean = mean(postList)
postSTD = std(postList)
postN = length(postList)

[h,p] = ttest(preList,postList)

Meandiff = postMean-preMean

