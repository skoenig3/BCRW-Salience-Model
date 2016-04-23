%% Code combine Statistcs across animals and image sets
load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\CombinedViewingBehavior.mat']);

tags = {'MP','TT','JN','IW'};
fixationdurations =[];
saccadedurations = [];
saccadeamplitudes = [];
numfixations = [];
distancebtwnfixations = [];
saccaderate = [];
for t= 1:length(tags);
    fixationdurations =[fixationdurations;allview{t}.fixduration];
    numfixations = [numfixations;sum(~isnan(allview{t}.fixduration),2)];
    saccadedurations = [saccadedurations;allview{t}.sacduration];
    saccadeamplitudes = [saccadeamplitudes;allview{t}.sacdist];
    distancebtwnfixations = [distancebtwnfixations;allview{t}.distbtwnfix];
    saccaderate = [saccaderate;allview{t}.timebtwfix];
end
%% plot saccade duration vs saccade arc length
saccadeamplitudes =  saccadeamplitudes(1:end)/24;
saccadeamplitudes(isnan( saccadeamplitudes)) = [];

saccadedurations = 5*saccadedurations(1:end);
saccadedurations(isnan(saccadedurations)) = [];

sthresh = prctile(saccadedurations,95);
sathresh = prctile(saccadeamplitudes,95);

rmd = find(saccadeamplitudes > sathresh);
saccadeamplitudes(rmd) = [];
saccadedurations(rmd) = [];
rmd = find(saccadedurations > sthresh);
saccadeamplitudes(rmd) = [];
saccadedurations(rmd) = [];

plot(saccadedurations,saccadeamplitudes,'.')
xlabel('Saccade Duration (ms)')
ylabel('Saccade Arc Length (dva)')
box off

xx = hist(saccadedurations,[10:5:90]);
xx = xx/max(xx)*max(saccadeamplitudes);
hold on
plot(10:5:90,xx,'r')
%%
for i = 10:5:90;
    ind = find(saccadedurations == i);
    plot(i,median(saccadeamplitudes(ind)),'k+','markersize',15);
end
%% durations in ms 5 ms per sample
fixationdurations = fixationdurations(1:end);
fixationdurations(isnan(fixationdurations)) = [];
fthresh = prctile(fixationdurations,95);
fixationdurations(fixationdurations > fthresh) = [];
meanfixationduration = mean(5*fixationdurations);
stdfixationdurations = std(5*fixationdurations);

saccadedurations = saccadedurations(1:end);
saccadedurations(isnan(saccadedurations)) = [];
sthresh = prctile(saccadedurations,95);
saccadedurations(saccadedurations > sthresh) = [];
meansaccadeduration = mean(5*saccadedurations);
stdsaccadedurations = std(5*saccadedurations);
%% number of fixations per image
meannumfixations = mean(numfixations);
stdnumfixations = std(numfixations);
%% saccade or fixaiton rate
saccaderate = saccaderate(1:end);
saccaderate(isnan(saccaderate)) = [];
srthresh = fthresh+sthresh;
saccaderate(saccaderate > srthresh) = [];
saccaderate = 1000./(5*saccaderate);
meansaccaderate = mean(saccaderate);
stdsaccaderate = std(saccaderate);
%% saccade amplitude and distance between fixatons in dva
saccadeamplitudes =  saccadeamplitudes(1:end);
saccadeamplitudes(isnan( saccadeamplitudes)) = [];
sathresh = prctile(saccadeamplitudes,95);
saccadeamplitudes(saccadeamplitudes > sathresh) = [];
meansaccadeamplitudes = mean(saccadeamplitudes)/24;
stdsaccadeamplitudes = std(saccadeamplitudes)/24;

distancebtwnfixations = distancebtwnfixations(1:end);
distancebtwnfixations(isnan(distancebtwnfixations)) = [];
meandistancebtwnfixations = mean(distancebtwnfixations)/24;
stddistancebtwnfixations = std(distancebtwnfixations)/24;
%%
clear,clc
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;

fixduration = NaN(2500,100);
sacduration = NaN(2500,100);
sacdist = NaN(2500,100);
numfixes = NaN(2500,1);
distbtwnfix = NaN(2500,100);
timebtwfix = NaN(2500,100);
count = 1;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    for mi = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{mi},'-MST-1.mat'); %must use MST-1.mat so won't pull MST-1.25 as well
        if ~isempty(str)
            for ii = 1:length(tags);
                strt = strfind(matfiles.mat{mi},tags{ii});
                if ~isempty(strt)
                    load(matfiles.mat{mi});
                    for cndlop = 1:2:length(fixationstats); %only uses novel viewing since images changes could alter natural behavior
                        reindexed = (cndlop+1)/2;
                        fixations = fixationstats{cndlop}.fixations;
                        fixationtimes = fixationstats{cndlop}.fixationtimes;
                        saccadetimes =  fixationstats{cndlop}.saccadetimes;
                        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                            fixations(:,1) = [];
                            fixationtimes(:,1) = [];
                        end
                        xy =fixationstats{cndlop}.XY;
                        sacduration(reindexed,1:length(saccadetimes)) = diff(saccadetimes,1)'+1;
                        fixduration(reindexed,1:length(fixationtimes)) = diff(fixationtimes,1)'+1;
                        numfixes(count) = size(fixationtimes,2);
                        for i = 1:size(fixations,2)-1
                            timebtwfix(count,i) = (fixationtimes(2,i+1)+fixationtimes(1,i+1))/2 ...
                                -(fixationtimes(2,i)+fixationtimes(1,i))/2;
                            if i == 1
                                x = fixations(1,i)-400;
                                y = fixations(2,i)-300;
                                distbtwnfix(count,i) = sqrt(x^2+y^2);
                                x = fixations(1,i+1)-fixations(1,i);
                                y = fixations(2,i+1)-fixations(2,i);
                                distbtwnfix(count,i+1) = sqrt(x^2+y^2);
                            else
                                x = fixations(1,i+1)-fixations(1,i);
                                y = fixations(2,i+1)-fixations(2,i);
                                distbtwnfix(count,i+1) = sqrt(x^2+y^2);
                            end
                        end
                        for i = 1:size(saccadetimes,2)
                            sacx = xy(1,saccadetimes(1,i):saccadetimes(2,i));
                            sacy = xy(2,saccadetimes(1,i):saccadetimes(2,i));
                            sacdist(count,i) = sum(sqrt(diff(sacx).^2+diff(sacy).^2));
                        end
                        count = count + 1;
                    end
                end
            end
        end
    end
end

%Fixation Duration
fixduration = fixduration(1:end);
fixduration(isnan(fixduration)) = [];
fthresh = prctile(fixduration,95);
fixduration(fixduration > fthresh) = [];
meanfixationduration = mean(5*fixduration);
stdfixationdurations = std(5*fixduration);
disp(['Fixation Duration ' num2str(meanfixationduration) ' ' setstr(177) ' ' num2str(stdfixationdurations)])
% Saccade Duration
saccadedurations = sacduration(1:end);
saccadedurations(isnan(saccadedurations)) = [];
sthresh =  prctile(saccadedurations,95);
saccadedurations(saccadedurations > sthresh) = [];
meansaccadeduration = mean(5*saccadedurations);
stdsaccadedurations = std(5*saccadedurations);
disp(['Saccade Duration ' num2str(meansaccadeduration) ' ' setstr(177) ' ' num2str(stdsaccadedurations)])
% saccade amplitude
saccadeamplitudes =  sacdist(1:end);
saccadeamplitudes(isnan( saccadeamplitudes)) = [];
sathresh = prctile(saccadeamplitudes,95);
saccadeamplitudes(saccadeamplitudes > sathresh) = [];
meansaccadeamplitudes = mean(saccadeamplitudes)/24;
stdsaccadeamplitudes = std(saccadeamplitudes)/24;
disp(['Saccade Arc Length ' num2str(meansaccadeamplitudes) ' ' setstr(177) ' ' num2str(stdsaccadeamplitudes)])
% Distance between fixations
distancebtwnfixations = distbtwnfix(1:end);
distancebtwnfixations(isnan(distancebtwnfixations)) = [];
meandistancebtwnfixations = mean(distancebtwnfixations)/24;
stddistancebtwnfixations = std(distancebtwnfixations)/24;
disp(['Distance btw fixations ' num2str(meandistancebtwnfixations) ' ' setstr(177) ' ' num2str(stddistancebtwnfixations)])
% number of fixations per image
meannumfixations = nanmean(numfixes);
stdnumfixations = nanstd(numfixes);
disp(['Number of fixations ' num2str(meannumfixations) ' ' setstr(177) ' ' num2str(stdnumfixations)])
% saccade or fixaiton rate
saccaderate = timebtwfix(1:end);
saccaderate(isnan(saccaderate)) = [];
srthresh = fthresh+sthresh;
saccaderate(saccaderate > srthresh) = [];
saccaderate = 1000./(5*saccaderate);
meansaccaderate = mean(saccaderate);
stdsaccaderate = std(saccaderate);
disp(['Saccade Rate ' num2str(meansaccaderate) ' ' setstr(177) ' ' num2str(stdsaccaderate)])