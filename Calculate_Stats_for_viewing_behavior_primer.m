%change location of CombinedViewing behavior data .m file
load('C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\CombinedViewingBehavior')
%% Get 1D PDFS
field = 'distbtwnfix';
name = 'Distance between Fixations (dva)';

bins = 100;
scale = 1/24;

cutoffl = 0;
cutoffu = 1000;
cutofflpct = 'microsaccades < 1 dva';
cutoffupct = 'unrealistic';

stats ={'','Mean','Median','STD','Min','Max','5%','25%','75%','95%'};

m1 = scale*getfield(allview{1},field);
m2 = scale*getfield(allview{2},field);
m3 = scale*getfield(allview{3},field);
m4 = scale*getfield(allview{4},field);
mall = [m1 m2 m3 m4];

m1(isnan(m1)) = [];
m2(isnan(m2)) = [];
m3(isnan(m3))= [];
m4(isnan(m4)) = [];
mall(isnan(mall)) = [];

means = [mean(m1) mean(m2) mean(m3) mean(m4) mean(mall)];
stds = [std(m1) std(m2) std(m3) std(m4) std(mall)];
medians = [median(m1) median(m2) median(m3) median(m4) median(mall)];
rngl = [min(m1) min(m2) min(m3) min(m4) min(mall)];
rngu = [max(m1) max(m2) max(m3) max(m4) max(mall)];

pct5 = [prctile(m1,5) prctile(m2,5) prctile(m3,5) prctile(m4,5) prctile(mall,5)];
pct25 = [prctile(m1,25) prctile(m2,25) prctile(m3,25) prctile(m4,25) prctile(mall,25)];
pct75 = [prctile(m1,75) prctile(m2,75) prctile(m3,75) prctile(m4,75) prctile(mall,75)];
pct95 = [prctile(m1,95) prctile(m2,95) prctile(m3,95) prctile(m4,95) prctile(mall,95)];

A = [means; medians; stds; rngl; rngu; pct5; pct25; pct75; pct95];
A = [stats' [[tags {'combined'}];num2cell(A)]];
xlswrite(field,A)

m1(m1 < cutoffl) = [];
m1(m1 > cutoffu) = [];
m2(m2 < cutoffl) = [];
m2(m2 > cutoffu) = [];
m3(m3 < cutoffl) = [];
m3(m3 > cutoffu) = [];
m4(m4 < cutoffl) = [];
m4(m4 > cutoffu) = [];
mall(mall < cutoffl) = [];
mall(mall > cutoffu) = [];

figure
hist(m1,bins);
title(tags{1})
xlabel(name);
box off

figure
hist(m2,bins);
title(tags{2})
xlabel(name);
box off

figure
hist(m3,bins);
title(tags{3})
xlabel(name);
box off

figure
hist(m4,bins);
title(tags{4})
xlabel(name);
box off

figure
hist(mall,bins);
title('Combined')
xlabel(name);
box off

means = [mean(m1) mean(m2) mean(m3) mean(m4) mean(mall)];
stds = [std(m1) std(m2) std(m3) std(m4) std(mall)];
medians = [median(m1) median(m2) median(m3) median(m4) median(mall)];
rngl = [min(m1) min(m2) min(m3) min(m4) min(mall)];
rngu = [max(m1) max(m2) max(m3) max(m4) max(mall)];

pct5 = [prctile(m1,5) prctile(m2,5) prctile(m3,5) prctile(m4,5) prctile(mall,5)];
pct25 = [prctile(m1,25) prctile(m2,25) prctile(m3,25) prctile(m4,25) prctile(mall,25)];
pct75 = [prctile(m1,75) prctile(m2,75) prctile(m3,75) prctile(m4,75) prctile(mall,75)];
pct95 = [prctile(m1,95) prctile(m2,95) prctile(m3,95) prctile(m4,95) prctile(mall,95)];

A = [means; medians; stds; rngl; rngu; pct5; pct25; pct75; pct95];
if cutofflpct ~= 0;
    if cutoffupct ~= 0
        A = [num2cell(A);cell(1,5); ...
            {'Lower cutoff:',cutoffl,'Upper Cutoff:',cutoffu,''};...
            {'Lower percentile:',cutofflpct,'Upper percentile:',cutoffupct,''}];
        extra = 1;
    else
        A = [num2cell(A);cell(1,5); ...
            {'Lower cutoff:',cutoffl,'','',''};...
            {'Lower percentile:',cutofflpct,'','',''}];
        extra = 1;
    end
elseif cutoffupct ~= 0;
    A = [num2cell(A);cell(1,5); ...
            {'','','Upper Cutoff:',cutoffu,''};...
            {'','','Upper percentile:',cutoffupct,''}];
        extra = 1;
else
    A = [num2cell(A);cell(1,5);{'Lower cutoff:',cutoffl,'Upper Cutoff:',cutoffu,''}];
    extra = 0;
end

A = [[stats'; cell(2+extra,1)] [[tags {'combined'}]; A]];
xlswrite([field '-filtered'],A)
%% Distance Profile
dm = {};
medianlen = 45;
for i = 1:4
    if size(allview{i}.distanceprofile,2) == medianlen;
        timewarp = 1:medianlen;
    else
        timewarp = round(linspace(1,size(allview{i}.distanceprofile,2),medianlen));
    end
    dm{i} = allview{i}.distanceprofile(:,timewarp);
end
dm{5} = [];
for i = 1:4
    dm{5} = [dm{5}; dm{i}];
end

for i = 1:4
figure
hold on
h = area(allview{i}.mediansac+1:size(dm{i},2),...
    max(nanmean(dm{i}/24+.1))*ones(1,size(dm{i},2)-allview{i}.mediansac));
set(h,'FaceColor',[.75 .75 .75])
set(h,'EdgeColor','none')
plot(nanmean(dm{i})/24,'linewidth',2);
hold off
title(tags{i})
xlabel('Warped Time')
ylabel('Distance (dva)')
ylim([min(nanmean(dm{i}/24)-0.1) max(nanmean(dm{i}/24+.1))])
end

figure
hold on
h = area(11:size(dm{5},2),...
    max(nanmean(dm{5}/24+.1))*ones(1,size(dm{5},2)-10));
set(h,'FaceColor',[.75 .75 .75])
set(h,'EdgeColor','none')
plot(nanmean(dm{5})/24,'linewidth',2);
hold off
title('Combined')
xlabel('Warped Time')
ylabel('Distance (dva)')
ylim([min(nanmean(dm{5}/24)-0.1) max(nanmean(dm{5}/24+.1))])
%% Angles
field =  'anglebtwfix';
tags = [tags {'Combined'}];

scale =  180/pi;%rads to degrees
n = (-180:180)*pi/180;

m{1} = scale*getfield(allview{1},field);
m{2} = scale*getfield(allview{2},field);
m{3} = scale*getfield(allview{3},field);
m{4} = scale*getfield(allview{4},field);
m{5} = [m{1} m{2} m{3} m{4}];

for i = 1:5
    [allprobsacangle] = hist(m{i}(~isnan(m{i})),360);
    allprobsacangle = [allprobsacangle(36:-1:1) allprobsacangle allprobsacangle(end:-1:end-36)];
    allprobsacangle = filtfilt(1/6*ones(1,6),1,allprobsacangle);
    allprobsacangle = allprobsacangle(37:end-37);
    allprobsacangle = allprobsacangle/sum(allprobsacangle);
    allprobsacangle = [allprobsacangle allprobsacangle(1)];
    
    figure
    polar(n,allprobsacangle)
    title(tags{i})
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    
end
%% 2D Fixation Location PDFs
f = fspecial('gaussian',[256,256],24);

m1 = allview{1}.densitymap;
m2 = allview{2}.densitymap;
m3 = allview{3}.densitymap;
m4 = allview{4}.densitymap;
mall = m1+m2+m3+m4;

m1 = imfilter(m1,f);
m2 = imfilter(m2,f);
m3 = imfilter(m3,f);
m4 = imfilter(m4,f);
mall = imfilter(mall,f);

figure
imagesc(m1)
axis off
title(tags{1})

figure
imagesc(m2)
axis off
title(tags{2})

figure
imagesc(m3)
axis off
title(tags{3})

figure
imagesc(m4)
axis off
title(tags{4})

figure
imagesc(mall)
axis off
title(tags{5})
%% Persistence
dm = {};
medianlen = 45;
for i = 1:4
    if size(allview{i}.persistence,2) == medianlen;
        timewarp = 1:medianlen;
    else
        timewarp = round(linspace(1,size(allview{i}.persistence,2),medianlen));
    end
    dm{i} = allview{i}.persistence(:,timewarp);
end
dm{5} = [];
for i = 1:4
    dm{5} = [dm{5}; dm{i}];
end

for i = 1:4
figure
hold on
h = area(allview{i}.mediansac+1:size(dm{i},2),...
    max(nanmean(dm{i}+.1))*ones(1,size(dm{i},2)-allview{i}.mediansac));
set(h,'FaceColor',[.75 .75 .75])
set(h,'EdgeColor','none')
plot(nanmean(dm{i}),'linewidth',2);
hold off
title(tags{i})
xlabel('Warped Time')
ylabel('Probability')
ylim([min(nanmean(dm{i})-0.1) max(nanmean(dm{i}+.1))])
end

figure
hold on
h = area(11:size(dm{5},2),...
    max(nanmean(dm{5}+.1))*ones(1,size(dm{5},2)-10));
set(h,'FaceColor',[.75 .75 .75])
set(h,'EdgeColor','none')
plot(nanmean(dm{5}),'linewidth',2);
hold off
title('Combined')
xlabel('Probability')
ylabel('Distance (dva)')
ylim([min(nanmean(dm{5})-0.1) max(nanmean(dm{5}+.1))])
%%
sacarclength = [];
sacamplitude = [];
for i = 1:4
    sacarclength = [sacarclength; allview{i}.sacdist/24];
    sacamplitude = [sacamplitude; allview{i}.sacamplitude/24];
end
sacarclength(isnan(sacarclength)) = [];
sacamplitude(isnan(sacamplitude)) = [];

sacamplitude(sacarclength > 42) = [];
sacarclength(sacarclength > 42) = [];

plot(sacamplitude,sacarclength,'.')
xlim([0 max(sacamplitude)])
ylim([0 42])
xlabel('Saccade Amplitude (dva)')
ylabel('Saccade Arc Length (dva)')

[r pvalue] = corrcoef(sacamplitude, sacarclength);
p = polyfit(sacamplitude, sacarclength,1);
%% Stats by Fixation and Saccade number
sacdist = [];
fixdur = [];
numfixes = [];
for i = 1:4
    sacdist = [sacdist; allview{i}.sacamplitude];
    fixdur = [fixdur; allview{i}.fixduration];
    temp = ~isnan(allview{i}.fixduration); 
    numfixes = [numfixes; reshape(sum(temp,2),[8,36])];
end
    
figure
errorbar(mean(numfixes),std(numfixes)/sqrt(8))%divide by number of image sets not # monkeys * # image sets
xlabel('Image Number')
ylabel('Number of Fixations')
box off

figure
errorbar(nanmean(sacdist/24),nanstd(sacdist/24)/sqrt(size(sacdist,1)))
xlabel('Saccade Number')
ylabel('Saccade Amplitude (dva)')
xlim([0 36])
box off

figure
errorbar(nanmean(fixdur*5),nanstd(fixdur*5)/sqrt(size(fixdur,1)))
xlabel('Fixation Number')
ylabel('Fixation Duration (ms)')
box off
xlim([0 36])
