function getViewingBehaviorCF(FIXATIONFILE,SAMPRATE,imageX,imageY,PLOTOPTIONS)
% Seth Koenig 06/06/2012 modified on 11/15/2012 for Fixations detected using
% Cluster Fixation (CF)

% Function extracts behavior from eye tracking data from cortex files for free
% viewing of natural scence. The function extracts % angles, distance, and 
% time between fixations.

%Inputs:
%   FIXATIONFILE:   Fixations extracted from cortex e.g. MP120606_1-fixation.mat
%   SAMPRATE:       Sampling rate of eye tracking data

%   PLOTOPTIONS: determines if function plots behavior statistics
%   PLOTOPTIONS = 'all' to plot all behavior statistics

%Outputs:
%   A file with the following name [FIXATIONFILE(1:end-13) '-ViewingBehavior'].
%   File contains the necessary variables and data for succesive functions
%   for modeling viewing behavior using a BCRW.
%   Saved variables are described by variable variablenames.

if nargin < 1
    error('Not enough inputs: function requires FIXATIONFILE');
end
if nargin < 2
    SAMPRATE = 5; %in ms/sampled point
end
if nargin < 3
    imageX = 800;
    imageY = 600;
end
if nargin < 5
    PLOTOPTIONS ='none';
end

transitionthreshold = 45; %counting as a turn/rotation
load(FIXATIONFILE);% loads fixaqtionstats

%-----Get fixations and saccade parameters----%
numtrials = round(length(fixationstats)/2);
fixduration = NaN(numtrials ,75);
distbtwnfix = NaN(numtrials ,75);
timebtwfix = NaN(numtrials,75);
sacduration = NaN(numtrials,75);
sacdist = NaN(numtrials,75);
sacamplitude = NaN(numtrials,75);
sacangle = NaN(numtrials,75);
sacangle_2fix = NaN(numtrials,75);
anglebtwfix = NaN(numtrials,75);
densitymap = zeros(600,800);
for cndlop = 1:2:length(fixationstats); %only uses novel viewing since images changes could alter natural behavior
    reindexed = (cndlop+1)/2;
    fixations = fixationstats{cndlop}.fixations;
    if ~isempty(fixations)
        fixationtimes = fixationstats{cndlop}.fixationtimes;
        saccadetimes =  fixationstats{cndlop}.saccadetimes;
        xy = fixationstats{cndlop}.XY;
        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
            fixations(:,1) = [];
            fixationtimes(:,1) = [];
        end
        xy =fixationstats{cndlop}.XY;
        sacduration(reindexed,1:length(saccadetimes)) = diff(saccadetimes,1)'+1;
        fixduration(reindexed,1:length(fixationtimes)) = diff(fixationtimes,1)'+1;
        for i = 1:size(fixations,2)-1
            fixx = round(fixations(1,i)); fixy = round(fixations(2,i));
            fixx(fixx < 1) = 1; fixx(fixx > 800) = 800;
            fixy(fixy < 1) = 1; fixy(fixy > 600) = 600;
            densitymap(fixy,fixx) = densitymap(fixy,fixx)+1;
            timebtwfix(reindexed,i) = (fixationtimes(2,i+1)+fixationtimes(1,i+1))/2 ...
                -(fixationtimes(2,i)+fixationtimes(1,i))/2; %time between middle of fixation and middle of next fixation
            if i == 1
                %distance from center crosshair @ 400,300
                x = fixations(1,i)-400; %change in x
                y = fixations(2,i)-300; %change in y
                distbtwnfix(reindexed,i) = sqrt(x^2+y^2);
                anglebtwfix(reindexed,i) = atan2(y,x);
                %distance from next fixation
                x = fixations(1,i+1)-fixations(1,i);
                y = fixations(2,i+1)-fixations(2,i);
                distbtwnfix(reindexed,i+1) = sqrt(x^2+y^2);
                anglebtwfix(reindexed,i+1) = atan2(y,x);
            else
                x = fixations(1,i+1)-fixations(1,i); %change in x
                y = fixations(2,i+1)-fixations(2,i); %change in y
                distbtwnfix(reindexed,i+1) = sqrt(x^2+y^2);
                anglebtwfix(reindexed,i+1) = atan2(y,x);
            end
        end
        for i = 1:size(saccadetimes,2)
            sacx = xy(1,saccadetimes(1,i):saccadetimes(2,i));
            sacy = xy(2,saccadetimes(1,i):saccadetimes(2,i));
            sacdist(reindexed,i) = sum(sqrt(diff(sacx).^2+diff(sacy).^2)); %arc length
            sacamplitude(reindexed,i) = sqrt((sacx(1)-sacx(end))^2+(sacy(1)-sacy(end))^2); %amplitude
            sacangle(reindexed,i) = atan2(diff(sacy(1:2)),diff(sacx(1:2))); %initial saccade angle
            sacangle_2fix(reindexed,i) =  atan2(diff(sacy(end-1:end)),diff(sacx(end-1:end))); %angle entering fixation
        end
    end
end

%---Angle between fixations---%
n = (-180:180)*pi/180;
[probanglebtwfix] = hist(anglebtwfix(~isnan(anglebtwfix)),360);
probanglebtwfix = [probanglebtwfix(36:-1:1) probanglebtwfix probanglebtwfix(end:-1:end-36)];
probanglebtwfix = filtfilt(1/18*ones(1,18),1,probanglebtwfix);
probanglebtwfix = probanglebtwfix(37:end-37);
probanglebtwfix = probanglebtwfix/sum(probanglebtwfix);
probanglebtwfix = [probanglebtwfix probanglebtwfix(1)];

%----Direction of Saccade leaving a fixation---%
[probsacangle] = hist(sacangle(~isnan(sacangle)),360);
probsacangle = [probsacangle(36:-1:1) probsacangle probsacangle(end:-1:end-36)];
probsacangle = filtfilt(1/18*ones(1,18),1,probsacangle);
probsacangle = probsacangle(37:end-37);
probsacangle = probsacangle/sum(probsacangle);
probsacangle = [probsacangle probsacangle(1)];

%---Parameter Profiles of Fixations and Saccades---%
variables = {'Dist','Vel','Accel','Rotation'};
fixlen = round(nanmedian(nanmedian(fixduration)))+10;
saclen = round(nanmedian(nanmedian(sacduration)))+10;
numfixes = sum(sum(~isnan(fixduration)));
numsacs = sum(sum(~isnan(sacduration)));
allfixations = zeros(numfixes,fixlen,length(variables));
allsaccades = zeros(numsacs,saclen,length(variables));
persistence.fix = zeros(numfixes,fixlen);
persistence.sac = zeros(numfixes,saclen);
distanceprofile.fix = NaN(numfixes,fixlen);
distanceprofile.sac = NaN(numfixes,saclen);
fixcount = 1;
saccount = 1;
for cndlop = 1:2:length(fixationstats); %only uses novel viewing since memory could alter natural behavior
    fixationtimes = fixationstats{cndlop}.fixationtimes;
    saccadetimes =  fixationstats{cndlop}.saccadetimes;
    xy =fixationstats{cndlop}.XY;
    if ~isempty(fixationtimes)
        if fixationtimes(1,1)<saccadetimes(1,1)
            fixcount2 = fixcount+1;
        else
            fixcount2 = fixcount;
        end
        for i = 1:size(saccadetimes,2);
            if i == 1
                x = xy(1,saccadetimes(1,i):saccadetimes(2,i)+7); %7 not 5 because going to index +2
                y = xy(2,saccadetimes(1,i):saccadetimes(2,i)+7); %7 not 5 because going to index +2
            elseif i == size(saccadetimes,2)
                x = xy(1,saccadetimes(1,i)-5:saccadetimes(2,i)); 
                y = xy(2,saccadetimes(1,i)-5:saccadetimes(2,i));
            else
                x = xy(1,saccadetimes(1,i)-5:saccadetimes(2,i)+7); %7 not 5 because going to index +2
                y = xy(2,saccadetimes(1,i)-5:saccadetimes(2,i)+7); %7 not 5 because going to index +2
            end 
            velx = diff(x);
            vely = diff(y);
            vel = sqrt(velx.^2+vely.^2);
            accel = abs(diff(vel));
            angle = 180*atan2(vely,velx)/pi;
            transitions = abs(diff(angle)) > transitionthreshold;
            vel = vel(1:end-1);
            rot = zeros(1,length(x)-2);
            dist = zeros(1,length(x)-2);
            for a = 1:length(x)-2;
                rot(a) = abs(angle(a)-angle(a+1));
                dist(a) = sqrt((x(a)-x(a+2)).^2 + (y(a)-y(a+2)).^2); %indexing +2 here
            end
            rot(rot > 180) = rot(rot > 180)-180;
            if  saclen == length(dist);
                timewarp = 1:length(dist);
            else
                timewarp = round(linspace(1,length(dist),saclen));
            end
            if i == size(saccadetimes,2)
                if fixationtimes(1,end) >= saccadetimes(2,end)
                    distanceprofile.sac(fixcount2,:) = vel(timewarp);
                    persistence.sac(fixcount2,:) = transitions(timewarp);
                end
            else
                distanceprofile.sac(fixcount2,:) = vel(timewarp);
                persistence.sac(fixcount2,:) = transitions(timewarp);
            end
            allsaccades(saccount,:,1) = dist(timewarp);
            allsaccades(saccount,:,2) = vel(timewarp);
            allsaccades(saccount,:,3) = accel(timewarp);
            allsaccades(saccount,:,4) = rot(timewarp);
            saccount = saccount+1;
            fixcount2 = fixcount2+1;
        end
        for i = 1:size(fixationtimes,2);
            if i == 1
                x = xy(1,fixationtimes(1,i):fixationtimes(2,i)+7);
                y = xy(2,fixationtimes(1,i):fixationtimes(2,i)+7);
            elseif i == size(fixationtimes,2)
                x = xy(1,fixationtimes(1,i)-5:fixationtimes(2,i));
                y = xy(2,fixationtimes(1,i)-5:fixationtimes(2,i));
            else
                x = xy(1,fixationtimes(1,i)-5:fixationtimes(2,i)+7);
                y = xy(2,fixationtimes(1,i)-5:fixationtimes(2,i)+7);
            end
            velx = diff(x);
            vely = diff(y);
            vel = sqrt(velx.^2+vely.^2);
            accel = abs(diff(vel));
            angle = 180*atan2(vely,velx)/pi;
            transitions = abs(diff(angle)) > transitionthreshold;
            vel = vel(1:end-1);
            rot = zeros(1,length(x)-2);
            dist = zeros(1,length(x)-2);
            for a = 1:length(x)-2;
                rot(a) = abs(angle(a)-angle(a+1));
                dist(a) = sqrt((x(a)-x(a+2)).^2 + (y(a)-y(a+2)).^2);
            end
            rot(rot > 180) = rot(rot > 180)-180;
            if  fixlen == length(dist);
                timewarp = 1:length(dist);
            else
                timewarp = round(linspace(1,length(dist),fixlen));
            end
            if i == 1
                if fixationtimes(1,1) >= saccadetimes(2,1)
                    distanceprofile.fix(fixcount,:) = vel(timewarp);
                    persistence.fix(fixcount,:) = transitions(timewarp);
                end
            else
                distanceprofile.fix(fixcount,:) = vel(timewarp);
                persistence.fix(fixcount,:) = transitions(timewarp);
            end
            allfixations(fixcount,:,1) = dist(timewarp);
            allfixations(fixcount,:,2) = vel(timewarp);
            allfixations(fixcount,:,3) = accel(timewarp);
            allfixations(fixcount,:,4) = rot(timewarp);
            fixcount = fixcount+1;
        end
    end
end
avgfixation= mean(allfixations,1);
avgfixprofile = zeros(size(avgfixation));
for i = 1:size(avgfixation,3);
    avgfixprofile(:,:,i) = filtfilt(1/3*ones(1,3),1,avgfixation(:,:,i));
    avgfixprofile(:,:,i) = avgfixprofile(:,:,i) - min(avgfixprofile(:,:,i));
    avgfixprofile(:,:,i) = avgfixprofile(:,:,i)/max(avgfixprofile(:,:,i));
end
avgsaccade= mean(allsaccades,1);
avgsacprofile = zeros(size(avgsaccade));
for i = 1:size(avgsaccade,3);
    avgsacprofile(:,:,i) = filtfilt(1/3*ones(1,3),1,avgsaccade(:,:,i));
    avgsacprofile(:,:,i) =  avgsacprofile(:,:,i) - min(avgsacprofile(:,:,i));
    avgsacprofile(:,:,i) = avgsacprofile(:,:,i)/max(avgsacprofile(:,:,i));
end

Statsbyfixation.fixatoinspertrial = sum(~isnan(fixduration),2)';
Statsbyfixation.meanfixationduration = SAMPRATE*nanmean(fixduration);
Statsbyfixation.stdfixationduration = SAMPRATE*nanstd(fixduration);
Statsbyfixation.numfix = sum(~isnan(fixduration));
Statsbyfixation.meansacdistance = nanmean(sacdist);
Statsbyfixation.stdsacdistance = nanstd(sacdist);
Statsbyfixation.numsacs = sum(~isnan(sacdist));

if strcmpi(PLOTOPTIONS,'all')
    figure
    hist(fixduration(~isnan(fixduration))*SAMPRATE,100)
    title('Fixation Duration')
    xlabel('Time (ms)')
    
    figure
    hist(distbtwnfix(~isnan(distbtwnfix)),100)
    title('Distance (Pixels)')
    
    figure
    hist(1000./(timebtwfix(~isnan(timebtwfix))*SAMPRATE),100)
    title('Fixation Rate')
    xlabel('Hz')
    
    figure
    hist(sacduration(~isnan(sacduration))*SAMPRATE,25)
    title('Saccade Duration')
    xlabel('Time (ms)')
    
    figure
    hist(sacdist(~isnan(sacdist)),100)
    title('Saccadic Distance')
    xlabel('Distance (Pixels)')
    
    figure
    subplot(1,2,1)
    polar(n,probsacangle)
    title('Probability Distribution of Saccade Angles Leaving a Fixation')
    
    subplot(1,2,2)
    polar(n,probanglebtwfix)
    title('Probability Distribution of Angles between Fixations')
    
    figure
    plot(sacdist(~isnan(sacdist)),sacduration(~isnan(sacduration))*SAMPRATE,...
        '*','markersize',3)
    title('Correlation between saccade duration and distance')
    xlabel('Distance (pixels)')
    ylabel('Duration (ms)')
    
    figure
    imagesc(densitymap)
    title('Probability Distribution of Fixations')
    axis off
    
    figure
    subplot(3,1,1)
    plot(Statsbyfixation.fixatoinspertrial);
    xlim([1 36])
    xlabel('Trial Number')
    ylabel('Number of Fixatons')
    title('Trials vs Number of Fixations')
    subplot(3,1,2)
    errorbar(Statsbyfixation.meanfixationduration,...
        Statsbyfixation.stdfixationduration./sqrt(Statsbyfixation.numfix))
    xl = find(Statsbyfixation.numfix < 5);
    xlim([0 xl(1)])
    xlabel('Fixation Number')
    ylabel('Fixation Duration (ms)')
    title('Fixation Duration vs Number of Fixations')
    subplot(3,1,3)
    errorbar(Statsbyfixation.meansacdistance,...
        Statsbyfixation.stdsacdistance./sqrt(Statsbyfixation.numsacs))
    xl = find(Statsbyfixation.numsacs < 5);
    xlim([0 xl(1)])
    xlabel('Saccade Number')
    ylabel('Saccade Distance (Pixels)')
    title('Saccade Distance vs Saccade Number')
    
    figure
    hold all
    h = area(5:fixlen-5,ones(1,fixlen-9));
    set(h,'FaceColor',[.75 .75 .75])
    set(h,'EdgeColor','none')
    for i =  1:size(avgfixprofile,3);
        plot(avgfixprofile(:,:,i),'linewidth',2)
    end
    hold off
    xlim([1 fixlen])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
    legend([{'fixation'} variables],'Location','NorthEastOutside');
    xlabel('Warped Time')
    ylabel('Normalized Value')
    title([FIXATIONFILE(1:10) '  Average-Smoothed Fixation Profile by Parameter'])
    
    figure
    hold all
    h(1) = area(1:5,ones(1,5));
    set(h(1),'FaceColor',[.75 .75 .75])
    set(h(1),'EdgeColor','none')
    h(2) = area(saclen-4:saclen,ones(1,5));
    set(h(2),'FaceColor',[.75 .75 .75])
    set(h(2),'EdgeColor','none')
    for i = 1:size(avgsacprofile,3)
        p(i) = plot(avgsacprofile(:,:,i),'linewidth',2);
    end
    hold off
    xlim([1 saclen])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
    legend([h(1) p],[{'fixation'} variables],'Location','NorthEastOutside');
    xlabel('Warped Time')
    ylabel('Normalized Value')
    title([FIXATIONFILE(1:10) '  Average-Smoothed Saccade Profile by Parameter'])
    
    figure
    hold on
    h(1) = area(1:5,ones(1,5));
    set(h(1),'FaceColor',[.75 .75 .75])
    set(h(1),'EdgeColor','none')
    h(2) = area(saclen-4:saclen,ones(1,5));
    set(h(2),'FaceColor',[.75 .75 .75])
    set(h(2),'EdgeColor','none')
    p = plot(sum(persistence.sac)./size(persistence.sac,1));
    hold off
    legend([h(1) p],{'fixation','persistence'},'Location','NorthEastOutside');
    xlabel('Warped Time')
    ylabel(['Probability of Saccade Angle Changeing > ' num2str(transitionthreshold) ' Degrees'])
end

variablenames={
    'fixduration: durations of fixations';
    'distbtwnfix: distance between fixations';
    'timebtwfix: time between fixations';
    'sacduration: saccadedurations';
    'sacdist: saccade are length';
    'sacamplitude: classic saccade amplitude';
    'sacangle: direction of saccade leaving a fixation';
    'sacangle_2fix: angle of saccade entering a fixation';
    'anglebtwfix: angles between successive fixations';
    'n: angles for probability distributions of angles';
    'probanglebtwfix: probability distribution of angles between successive fixations';
    'probsacangle: probability distrubution of direction of saccade leaving a fixation';
    'persistence: persistence of eye moving in previous direction (change < 45 degrees)';
    'densitymap: probability distribtuion for locations of fixations';
    'variables: variables for fixation and saccade profiles';
    'allfixations: all variables by fixation';
    'allsaccades: all variables by saccade';
    'avgfixprofile: average variable profile of all fixations warped to median fixation duration';
    'avgsacprofile: average variable profile of all saccades warped to median saccade duration';
    'distanceprofile: profile of distances by paired saccaed and fixation';
    'Statsbyfixation: stats for num fixation per trial, fixation duration by fixation#, & sac distance by saccade #';
    };

save([FIXATIONFILE(1:end-13) '-ViewingBehavior'],'fixduration',...
    'distbtwnfix','timebtwfix','sacduration','sacdist','sacangle','anglebtwfix',...
    'densitymap','n','probanglebtwfix','probsacangle','persistence','variables',...
    'allfixations','allsaccades','avgfixprofile','avgsacprofile',...
    'Statsbyfixation','distanceprofile','sacangle_2fix','sacamplitude')