%% Code Below Runs the core portion of the BCRW-Salience Model %%
% Seth Koenig 2012-2015. skoenig@gatech.edu

% This code does the following...search for each section using "[#]"
% 1. Create Salience Maps ~ 2 hours for 8 image sets of 36 image each
% 2. Extracts Fixations and Saccade times using Cluster Fix
% 3. Calculates the Salience Value at Each Fixation
% 4. Calculate Average Salience at Each Fixation Across mutliple data sets
% 5. Extract Viewing Behavior Statistics by image set and monkey
% 6. Combine Viewing Behavior by Monkey
% 7. Determine if Salence alters IOR (inhibition of return) or vice a versa
% 8. Run BCRW simulations ~24-36 hours for 8 image sets
% 9. Calculate goodness of fit of BCRW for fixation location using KL-Divergence
% 10. Calculate goodness of fit of BCRW for fixation location using AUC ROC
% 11. Calculate goodness of fit of BCRW for fixation order using AUC ROC

% REQUIRED functions-Core functions do most of the math and simulation
%   1) bin2: creates binned matrices
%   2) ClusterFixation: detects fixations and saccades using k-means cluster 
%   analysis. It is 100% essential to use this algorithm or an algorithm that
%   ouputs data in a similar structure.
%   3) fixationSalience_and_significanceCF: code determines salience,
%   salience contrast, and image intensity values at fixation locations.
%   4) getSalienceMap: creates saliency maps based on Itti et al 1998.
%       a) can use  getSalienceMapgray to create saliency maps for RGB images
%       that are essential gray scale.
%   5) getViewingBehvaiorCF: calculates viewing behavior statitics
%   6) runBCRW_CF: code runs the BCRW
%   7) SalienceIOR: calculates rate, duration, and salience
%   at return fixations and compares.


% RECOMMENDED functions/code:
%   1) AverageMaps: calculates the average salience, BCRW, and image
%   intensity maps
%   Calculates goodness of fit using KL divergence.
%   2) CombinedSaccadeandFixationProfiles: calculates the average fixation
%   and saccade parameter (distnace, velocity, acceleraiton, and rotaion)
%   values from data found by ClusterFixation.
%   3) FixationDetectionUsingVel_Accel_Thresholds: comparable fixation
%   detection algorithm to Cluster_Fixation but using velocity and
%   accelration thresholds.
%   4) fixationSalience_and_significanceCF_contrast: calculate the values
%   of each cotrast type that contributes to the saliency map at fixation
%   locations.
%   5) run_BCRW_parameter_sweep: runs a parameter sweep to determine the value
%   of 4 parameters not derived from actual viewing behaivor statistics.
%   6) BCRW_paramter_sweep: calls run_BCRW_parameter_sweep and analyzes the
%   data

% IMPORTANT NOTES: This code is totally automated. You may however need to edit
% lines of code such as image_sets and scm_image_dir so that you can analyze
% your data as it is organized on your computer. I suggest organizing your data
% in the following way. For each SCM (scene manipulation) image set create
% a new folder and place the images for that SCM set and the cortex data
% files for every monkey (e.g 'MP100807.2') that corresponds to that image
% folder. Create these folders for each SCM set--the 'image_sets' variable
% contains the names for these SCM set folders. The location of all these
% image sets is to be saved under the 'scm_image_dir' variable. When for
% example the saliency map code runs the a matlab file containing the saliency
% map will be placed in the SCM folder where the image was found. Also, make
% sure the save command saves to the desired folder. I suggest using the search
% function to fix these. Additionally, many of the sub-functions assume the
% structure of SCM images where novel images are followed by the
% presentation of the same image, the image with a replaced object, or the
% image with a moved object. Please be aware the code will have to be
% changed to use data organized in a different fashion or analyze data from
% familar images. FINALLY THIS CODE USES INVERTED IMAGE INTENSITIES in
% other words based on statistical results the code assumes lower image
% intensities predict fixations better than high image intensities.

%---[1] Get the salience maps for multiple task sets---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
for imset = 1:length(image_sets);
    dirName = [scm_image_dir image_sets{imset}];
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    imageindexes = [];
    for i = 1:length(fileList)
        bmps = strfind(fileList{i},'bmp');
        if ~isempty(bmps)
            if double(fileList{i}(bmps-2)) <= 57 %ascii for number
                imageindexes = [imageindexes i];
            end
        end
    end
    for i = 1:length(imageindexes)
        imagefile = fileList{imageindexes(i)};
        getSalienceMap(imagefile)
    end
end
%%
%---[2] Get Fixations and Saccades from Behavioral Files---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
for imset = 1:length(image_sets);
    dirName = [scm_image_dir image_sets{imset}];
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    eyeindexes = [];
    for i = 1:length(fileList)
        period = strfind(fileList{i},'.');
        if (length(fileList{i}) == period(end)+1) && ...
                (double(fileList{i}(period+1)) <=57)
            eyeindexes = [eyeindexes i];
        end
    end
    for i = 1:length(eyeindexes)
        cortexfile = fileList{eyeindexes(i)};
        [fixationstats,trialtype] = ClusterFixation(cortexfile);
        save([cortexfile(1:end-2) '_' cortexfile(end) '-fixation.mat'],'fixationstats','trialtype')
    end
end
%%
%---[3] Calculate Salience at Each Fixation---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
PLOTOPTIONS = 'none';
imageX = 800; imageY = 600;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    for eyefile = eyedatafiles;
        fixationSalience_and_significanceCF(matfiles.mat{eyefile},imageX,imageY,PLOTOPTIONS)
    end
end
%%
%%---[4] Calculate Average Salience at Each Fixation Across mutliple data sets---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
minlen = 100;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'FixationStatistics.mat');
        if ~isempty(str)
            statfiles = [statfiles i];
        end
    end
    for stat = statfiles;
        load(matfiles.mat{stat},'statistics')
        minlen = min(minlen,size(statistics.numbervalues,3));
    end
end

%override minlen above. Code wasn't written to get the median number of
%fixations. Minlen works well for plotting though by itself. Was 34. SDK
minlen = 40; %median number of fixations according to viewing behavior section. 

alldata = NaN(36*length(image_sets)*4,3,2,minlen);
allshuffled = NaN(36*length(image_sets)*4,3,2,minlen);
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    statfiles = zeros(1,length(tags));
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'FixationStatistics.mat'));
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    statfiles(ii) = i;
                end
            end
        end
    end
    
    for stat = 1:4; %iterate through monkey
        if statfiles(stat) ~= 0
            i = length(tags)*(SET-1)+stat;
            load(matfiles.mat{statfiles(stat)})
            combineddata = shuffunshuffdata{2}{3};
            combinedshuffled = shuffunshuffdata{1}{3};
            alldata(36*(i-1)+1:i*36,:,:,:) = combineddata(:,:,:,1:minlen);
            allshuffled(36*(i-1)+1:i*36,:,:,:) = combinedshuffled(:,:,:,1:minlen);
        end
    end
end

allmeanvals = NaN(3,2,minlen);
allstdvals = NaN(3,2,minlen);
allnumvals = NaN(3,2,minlen);
for i = 1:size(allmeanvals,1)
    for ii = 1:size(allmeanvals,2)
        for iii = 1:minlen
            allmeanvals(i,ii,iii) = nanmean(alldata(:,i,ii,iii));
            allstdvals(i,ii,iii) = nanstd(alldata(:,i,ii,iii));
            allnumvals(i,ii,iii) = sum(~isnan(alldata(:,i,ii,iii)));
        end
    end
end

% z-test of means agains random distributions assuming mean is larger
allzp = NaN(size(allmeanvals)); %p-values
allcI = NaN(size(allmeanvals)); %top confidence interval value, lowest is typticall 0/-Inf
for i = 1:size(allmeanvals,1)
    for ii = 1:size(allmeanvals,2)
        shuffledvals = allshuffled(:,i,ii,:);
        shuffledvals(isnan(shuffledvals)) = [];
        for iii = 1:size(allmeanvals,3)
            [~,p,ci] = ztest(shuffledvals,allmeanvals(i,ii,iii),std(shuffledvals),...
                0.05);
            allzp(i,ii,iii) = p;
            if i == 3;
                allcI(i,ii,iii) = ci(1);
            else
                allcI(i,ii,iii) = ci(2);
            end
        end
    end
end

allstatistics.meanvalues = allmeanvals;
allstatistics.stdvalues = allstdvals;
allstatistics.numbervalues = allnumvals;
allstatistics.pvalues = allzp;
allstatistics.confidenceintervals = allcI;

clrs = ['rb'];
legendlabels = {'Salience','Salience Contrast','Image Intensity'};
averaginglabels = {' at mean fixation location',' mean during fixation'};
for i = 1:size(allcI,1)
    figure
    hold on
    for ii = 1:size(allcI,2)
        p(ii)= plot(allcI(i,ii)*ones(1,size(allcI,3)),['--' clrs(ii)]);
        pp(ii) = errorbar(shiftdim(allmeanvals(i,ii,:)),...
            shiftdim(allstdvals(i,ii,:))./sqrt(shiftdim(allnumvals(i,ii,:))),clrs(ii));
    end
    hold off
    title([legendlabels{i} ' at a fixation by fixation number across all'...
        ' all monkeys and image sets'])
    legend([p pp],{['Chance ' legendlabels{i} averaginglabels{1}],...
        ['Chance ' legendlabels{i} averaginglabels{2}],...
        [legendlabels{i} averaginglabels{1}], ...
        [legendlabels{i} averaginglabels{2}]},'Location','NorthEastOutside')
    xlabel('Fixation Number')
    ylabel('Normalized Value')
end

allstatvariablenames = {
    'alldata: cell arrray containing combined values for salience, salience';
    'contrast, and image intensity at fixations across all monkeys and data files.';
    'Each of these cells are arranged by row,column,and z-column in the following';
    'manner:  Rows are arranged as Salience, salience contrast, and intensity;';
    'Columns indicate if data is these parmaters at the average fixation cooridante';
    '(col 1) or the average of the parameters during a fixatoin; Z-column is';
    'organized by fixation number';
    '';
    'allshuffled: same as alldata but shuffled data';
    '';
    'allstatistics: structure with results from a z-test to determine if fixations';
    'occur at salience, salience contrasts, and image intesntisy valeus at rates';...
    'higher than what would be expected by chance. Random distributions from each';
    'parameter is compared to the mean value at that parameter by fixation';
    'statistics contains means, std, number of fixations, p-values,and confidence intervals';
    'These variables are arranged by row,column,z-column. Rows are arranged as';
    'Salience, salience contrast, and intensity. Columns indicate if data is these';
    'paramaters at the average fixation cooridante (col 1) or the average of the';
    'parameters during a fixation. Z-column is organized by fixation number';
    };

save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\CombinedSalienceStatistics.mat'],'alldata','allshuffled',...
    'allstatistics','allstatvariablenames');
%%
%---[5] Extract Viewing Behavior---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
PLOTOPTIONS = 'none';
imageX = 800;
imageY = 600;
SAMPRATE = 5;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    for eyefile = eyedatafiles;
        getViewingBehaviorCF(matfiles.mat{eyefile},SAMPRATE,imageX,imageY,PLOTOPTIONS)
    end
end
%%
%---[6] Combine Viewing Behavior by Monkey---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
medianfix = NaN(length(image_sets),length(tags));
mediansac = NaN(length(image_sets),length(tags));
mediannumfix = NaN(length(image_sets),length(tags));
mediannumsac = NaN(length(image_sets),length(tags));
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'ViewingBehavior');
        if ~isempty(str)
            for ii = 1:length(tags);
                strt = strfind(matfiles.mat{i},tags{ii});
                if ~isempty(strt)
                    load(matfiles.mat{i},'fixduration','sacamplitude','avgsacprofile','avgfixprofile');
                    mediannumfix(SET,ii) = median(sum(~isnan(fixduration')));
                    mediannumsac(SET,ii) = median(sum(~isnan(sacamplitude')));
                    medianfix(SET,ii) = size(avgfixprofile,2);
                    mediansac(SET,ii) = size(avgsacprofile,2);
                end
            end
        end
    end
end
%%
%get the median number of fixations and saccades across all monkeys and sets
mediannumfix = round(median(mediannumfix)); 
mediannumsac = round(median(mediannumsac));
%get the median "duration" of time warped data
medianfix = round(median(medianfix));
mediansac = round(median(mediansac));

allview = cell(1,length(tags));
for i = 1:length(tags)
    allview{i}.densitymap = zeros(600,800);
    allview{i}.allfixations = [];
    allview{i}.allsaccades = [];
    allview{i}.persistence =[];
    allview{i}.anglebtwfix = [];
    allview{i}.sacangle_2fix = [];
    allview{i}.distanceprofile = [];
    allview{i}.distbtwnfix = [];
    allview{i}.fixduration = [];
    allview{i}.sacangle = [];
    allview{i}.sacdist = [];
    allview{i}.sacamplitude = [];
    allview{i}.sacduration = [];
    allview{i}.timebtwfix = [];
end

for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'ViewingBehavior');
        if ~isempty(str)
            for ii = 1:length(tags);
                strt = strfind(matfiles.mat{i},tags{ii});
                if ~isempty(strt)
                    load(matfiles.mat{i});
                    if size(allfixations,2) == medianfix(ii);
                        timewarp = 1:size(allfixations,2);
                    else
                        timewarp = round(linspace(1,size(allfixations,2),medianfix(ii)));
                    end
                    distanceprofile.fix = distanceprofile.fix(:,timewarp);
                    persistence.fix = persistence.fix(:,timewarp);
                    persistence.fix = persistence.fix(:,6:end-5); %first and last 5 are extended past eye movement
                    distanceprofile.fix = distanceprofile.fix(:,6:end-5);  %first and last 5 are extended past eye movement
                    allview{ii}.allfixations = [allview{ii}.allfixations;...
                        allfixations(:,timewarp,:)];
                    if size(allsaccades,2) == mediansac(ii);
                        timewarp = 1:size(allsaccades,2);
                    else
                        timewarp = round(linspace(1,size(allsaccades,2),mediansac(ii)));
                    end
                    allview{ii}.allsaccades = [allview{ii}.allsaccades;...
                        allsaccades(:,timewarp,:)];
                    distanceprofile.sac = distanceprofile.sac(:,timewarp);
                    distanceprofile.sac = distanceprofile.sac(:,6:end-5);  %first and last 5 are extended past eye movement
                    persistence.sac = persistence.sac(:,timewarp);
                    persistence.sac = persistence.sac(:,6:end-5); %first and last 5 are extended past eye movement
                    allview{ii}.persistence = [ allview{ii}.persistence;
                        [persistence.sac persistence.fix]];
                    allview{ii}.anglebtwfix = [allview{ii}.anglebtwfix;anglebtwfix];
                    allview{ii}.sacangle_2fix = [allview{ii}.sacangle_2fix;...
                        sacangle_2fix];
                    allview{ii}.densitymap = allview{ii}.densitymap+densitymap;
                    allview{ii}.distanceprofile = [allview{ii}.distanceprofile;
                        [distanceprofile.sac distanceprofile.fix]];
                    allview{ii}.distbtwnfix = [allview{ii}.distbtwnfix;distbtwnfix];
                    allview{ii}.fixduration = [allview{ii}.fixduration;...
                        fixduration];
                    allview{ii}.sacangle = [allview{ii}.sacangle;sacangle];
                    allview{ii}.sacdist = [allview{ii}.sacdist;sacdist];
                    allview{ii}.sacamplitude = [allview{ii}.sacamplitude;sacamplitude];
                    allview{ii}.sacduration = [allview{ii}.sacduration;sacduration];
                    allview{ii}.timebtwfix = [allview{ii}.timebtwfix;...
                        timebtwfix];
                    allview{ii}.mediansac = mediansac(ii)-10; %reduce since the 1st 5 and last 5 are data from the other eye movement
                    allview{ii}.medianfix = medianfix(ii)-10; %reduce since the 1st 5 and last 5 are data from the other eye movement
                    allview{ii}.mediannumfix = mediannumfix(ii);
                    allview{ii}.mediannumsac = mediannumsac(ii);
                end
            end
        end
    end
end

clearvars -except scm_image_dir image_sets allview tags

SAMPRATE = 5;
n = (-180:180)*pi/180;
variables = {'Dist','vel','accel','rot'};
f = fspecial('gaussian',[256,256],24);
graphnum = gcf;
if graphnum == 1
    graphnum = 0;
end

for i = 1:length(tags)
    
    
        %---Stats by fixation---%
    allStatsbyfixation{i}.fixatoinspertrial = reshape(sum(~isnan(allview{i}.fixduration),2),[],36);
    allStatsbyfixation{i}.meanfixationduration = SAMPRATE*nanmean(allview{i}.fixduration);
    allStatsbyfixation{i}.stdfixationduration = SAMPRATE*nanstd(allview{i}.fixduration);
    allStatsbyfixation{i}.numfix = sum(~isnan(allview{i}.fixduration));
    allStatsbyfixation{i}.meansacdistance = nanmean(allview{i}.sacdist);
    allStatsbyfixation{i}.stdsacdistance = nanstd(allview{i}.sacdist);
    allStatsbyfixation{i}.numsacs = sum(~isnan(allview{i}.sacdist));
    
    [allprobanglebtwfix] = hist(allview{i}.anglebtwfix(~isnan(allview{i}.anglebtwfix)),360);
    allprobanglebtwfix = [allprobanglebtwfix(36:-1:1) allprobanglebtwfix allprobanglebtwfix(end:-1:end-36)];
    allprobanglebtwfix = filtfilt(1/6*ones(1,6),1,allprobanglebtwfix);
    allprobanglebtwfix = allprobanglebtwfix(37:end-37);
    allprobanglebtwfix = allprobanglebtwfix/sum(allprobanglebtwfix);
    allprobanglebtwfix = [allprobanglebtwfix allprobanglebtwfix(1)];
    
    figure(graphnum+1)
    subtitle('Distribution of angles between fixations')
    hax(1,i) = subplot(2,2,i);
    polar(n,allprobanglebtwfix)
    title(tags{i})
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    
    [allprobsacangle] = hist(allview{i}.sacangle(~isnan(allview{i}.sacangle)),360);
    allprobsacangle = [allprobsacangle(36:-1:1) allprobsacangle allprobsacangle(end:-1:end-36)];
    allprobsacangle = filtfilt(1/6*ones(1,6),1,allprobsacangle);
    allprobsacangle = allprobsacangle(37:end-37);
    allprobsacangle = allprobsacangle/sum(allprobsacangle);
    allprobsacangle = [allprobsacangle allprobsacangle(1)];
    
    figure(graphnum+2)
    subtitle('Distribution of angles leaving a fixation')
    hax(2,i) = subplot(2,2,i);
    polar(n,allprobsacangle)
    title(tags{i})
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    
    %---Stats by fixation---%
    allStatsbyfixation{i}.fixatoinspertrial = reshape(sum(~isnan(allview{i}.fixduration),2),[],36);
    allStatsbyfixation{i}.meanfixationduration = SAMPRATE*nanmean(allview{i}.fixduration);
    allStatsbyfixation{i}.stdfixationduration = SAMPRATE*nanstd(allview{i}.fixduration);
    allStatsbyfixation{i}.numfix = sum(~isnan(allview{i}.fixduration));
    allStatsbyfixation{i}.meansacdistance = nanmean(allview{i}.sacdist);
    allStatsbyfixation{i}.stdsacdistance = nanstd(allview{i}.sacdist);
    allStatsbyfixation{i}.numsacs = sum(~isnan(allview{i}.sacdist));
    
    figure(graphnum+3)
    subtitle('Distribution of Fixation Durations')
    hax(3,i) = subplot(2,2,i);
    fixduration = allview{i}.fixduration(~isnan(allview{i}.fixduration))*SAMPRATE;
    fixduration(fixduration > 500) = [];
    hist(fixduration,95)
    xlabel('Time (ms)')
    title(tags{i})
    
    figure(graphnum+4)
    subtitle('Distribution of Distances between Fixations')
    hax(4,i) = subplot(2,2,i);
    hist(allview{i}.distbtwnfix(~isnan(allview{i}.distbtwnfix)),100)
    xlabel('Distance (Pixels)')
    title(tags{i})
    
    figure(graphnum+5)
    subtitle('Fixation Rate')
    hax(5,i) = subplot(2,2,i);
    hist(1000./(allview{i}.timebtwfix(~isnan(allview{i}.timebtwfix))*SAMPRATE),100)
    xlabel('Hz')
    title(tags{i})
    
    figure(graphnum+6)
    subtitle('Distribution of Saccade Durations (Arc Lenght)')
    hax(6,i) = subplot(2,2,i);
    sacduration = allview{i}.sacduration(~isnan(allview{i}.sacduration))*SAMPRATE;
    sacduration(sacduration > 150) = [];
    hist(sacduration,28)
    xlim([0 100])
    xlabel('Time (ms)')
    title(tags{i})
    
    
    figure(graphnum+7)
    subtitle('Distribution of Saccade Distances')
    hax(7,i) = subplot(2,2,i);
    sacdist = allview{i}.sacdist(~isnan(allview{i}.sacdist));
    sacdist(sacdist > 800) = [];
    hist(sacdist,100)
    xlabel('Distance (Pixels)')
    title(tags{i})
    
    figure(graphnum+8)
    subtitle('Correlation between saccade duration and distance')
    hax(8,i) = subplot(2,2,i);
    plot(allview{i}.sacdist(~isnan(allview{i}.sacdist)),...
        allview{i}.sacduration(~isnan(allview{i}.sacduration))*SAMPRATE,...
        '*','markersize',2)
    title(tags{i})
    xlabel('Distance (pixels)')
    ylabel('Duration (ms)')
    xlim([0 1000])
    ylim([0 200])
    
    figure(graphnum+9)
    subtitle('Probability Distribution of Fixations')
    hax(9,i) = subplot(2,2,i);
    densitymap = allview{i}.densitymap;
    densitymap = imfilter(densitymap,f);
    densitymap = densitymap./sum(sum(densitymap));
    imagesc(densitymap)
    title(tags{i})
    axis off
    
    figure(graphnum+10)
    subtitle('Fixation Statistics by Trial Number, Fixation Number, or Saccade Number')
    hax(10,i) = subplot(2,2,i);
    title(tags{i})
    hold all
    errorbar(mean(allStatsbyfixation{i}.fixatoinspertrial),...
        std(allStatsbyfixation{i}.fixatoinspertrial)/...
        size(allStatsbyfixation{i}.fixatoinspertrial,1));
    errorbar(allStatsbyfixation{i}.meanfixationduration,...
        allStatsbyfixation{i}.stdfixationduration./sqrt(allStatsbyfixation{i}.numfix))
    xl = find(allStatsbyfixation{i}.numfix < 50);
    xlim([0 xl(1)])
    ylim([0 400])
    errorbar(allStatsbyfixation{i}.meansacdistance,...
        allStatsbyfixation{i}.stdsacdistance./sqrt(allStatsbyfixation{i}.numsacs))
    if i == 1;
        ylabel('Number of Fixations, Fixation Duration (ms),Saccade Distance (Pixels)')
        set(get(gca,'YLabel'),'Position',[-7 -0.3 0])
    end
    if i == 2
        legend('# Fixations','Fixation Duration','Saccade Distance','Location',...
            'NorthEastOutside')
    end
    if i > 2
        xlabel('Trial Number, Fixation Number, or Saccade Number')
    end
    
    avgfixation= mean(allview{i}.allfixations,1);
    fixlen = size(avgfixation,2);
    avgfixprofile = zeros(size(avgfixation));
    for ii = 1:size(avgfixation,3);
        avgfixprofile(:,:,ii) = filtfilt(1/3*ones(1,3),1,avgfixation(:,:,ii));
        avgfixprofile(:,:,ii) = avgfixprofile(:,:,ii) - min(avgfixprofile(:,:,ii));
        avgfixprofile(:,:,ii) = avgfixprofile(:,:,ii)/max(avgfixprofile(:,:,ii));
    end
    avgsaccade= mean(allview{i}.allsaccades,1);
    saclen = size(avgsaccade,2);
    avgsacprofile = zeros(size(avgsaccade));
    for ii = 1:size(avgsaccade,3);
        avgsacprofile(:,:,ii) = filtfilt(1/3*ones(1,3),1,avgsaccade(:,:,ii));
        avgsacprofile(:,:,ii) =  avgsacprofile(:,:,ii) - min(avgsacprofile(:,:,ii));
        avgsacprofile(:,:,ii) = avgsacprofile(:,:,ii)/max(avgsacprofile(:,:,ii));
    end
    
    figure(graphnum+11)
    subtitle('Average-Smoothed Fixation Profile by Parameter')
    hax(11,i) = subplot(2,2,i);
    title(tags{i})
    hold all
    h = area(5:fixlen-5,ones(1,fixlen-9));
    set(h,'FaceColor',[.75 .75 .75])
    set(h,'EdgeColor','none')
    for ii =  1:size(avgfixprofile,3);
        plot(avgfixprofile(:,:,ii),'linewidth',2)
    end
    hold off
    xlim([1 fixlen])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
    if i == 1
        ylabel('Normalized Value')
        set(get(gca,'YLabel'),'Position',[-2 -0.2 0])
    end
    if i == 2
        legend([{'fixation'} variables],'Location','NorthEastOutside');
    end
    if i > 2
        xlabel('Warped Time')
    end
    
    figure(graphnum+12)
    subtitle('Average-Smoothed Saccade Profile by Parameter')
    hax(12,i) = subplot(2,2,i);
    title(tags{i})
    hold all
    h1 = area(1:5,ones(1,5));
    set(h1,'FaceColor',[.75 .75 .75])
    set(h1,'EdgeColor','none')
    h2 = area(saclen-4:saclen,ones(1,5));
    set(h2,'FaceColor',[.75 .75 .75])
    set(h2,'EdgeColor','none')
    for ii = 1:size(avgsacprofile,3)
        p(ii) = plot(avgsacprofile(:,:,ii),'linewidth',2);
    end
    hold off
    xlim([1 saclen])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
    if i == 1
        ylabel('Normalized Value')
        set(get(gca,'YLabel'),'Position',[-2 -0.2 0])
    end
    if i == 2
        legend([h1 p],[{'fixation'} variables],'Location','NorthEastOutside');
    end
    if i > 2
        xlabel('Warped Time')
    end
    
    figure(graphnum+13)
    subtitle('Probability of Saccade Angle Changeing > 45 Degrees')
    hax(13,i) = subplot(2,2,i);
    title(tags{i})
    hold on
    h  = area(allview{i}.mediansac+1:size(allview{i}.persistence,2),...
        ones(1,size(allview{i}.persistence,2)-allview{i}.mediansac));
    set(h,'FaceColor',[.75 .75 .75])
    set(h,'EdgeColor','none')
    p = plot(mean(allview{i}.persistence));
    hold off
    xlim([1 size(allview{i}.persistence,2)])
    if i == 1
        ylabel('Probability of Saccade Angle Changeing > 45 Degrees')
        set(get(gca,'YLabel'),'Position',[-3 -0.3 0])
    end
    if i == 2
        legend([h p],{'fixation','persistence'},'Location','NorthEastOutside');
    end
    if i > 2
        xlabel('Warped Time')
    end
    
    figure(graphnum+14)
    subtitle('Saccade and Fixation Distance')
    hax(14,i) = subplot(2,2,i);
    plot(nanmean(allview{i}.distanceprofile))
    title(tags{i})
    xlabel('Warped Time')
    ylabel('Distance (pixels)')
    
    [allprobangle2fix] = hist(allview{i}.sacangle_2fix(~isnan(allview{i}.sacangle_2fix)),360);
    allprobangle2fix = [allprobangle2fix(36:-1:1) allprobangle2fix allprobangle2fix(end:-1:end-36)];
    allprobangle2fix = filtfilt(1/6*ones(1,6),1,allprobangle2fix);
    allprobangle2fix = allprobangle2fix(37:end-37);
    allprobangle2fix = allprobangle2fix/sum(allprobangle2fix);
    allprobangle2fix = [allprobangle2fix allprobangle2fix(1)];
    
    figure(graphnum+15)
    subtitle('Distribution of saccade angles entering fixations')
    hax(15,i) = subplot(2,2,i);
    polar(n,allprobangle2fix)
    title(tags{i})
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    
    figure(graphnum+16)
    subtitle('Distribution of Saccade Amplitudes')
    hax(16,i) = subplot(2,2,i);
    sacdist = allview{i}.sacamplitude(~isnan(allview{i}.sacamplitude));
    sacdist(sacdist > 800) = [];
    hist(sacdist,100)
    xlabel('Distance (Pixels)')
    title(tags{i})
end

screen_size = get(0, 'ScreenSize');
figdir = ['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'Plots-All SCM Files\Viewing Behavior\'];
figuretitles = {
    'Distribution of Angles between fixations';
    'Distribution of Angles Leaving a fixation';
    'Distribution of Fixation Durations';
    'Distribution of Distances Between Fixations';
    'Fixation (Saccade) Rate';
    'Distribution of Saccade Durations';
    'Distribution of Saccade Distances';
    'Correlation between Saccade Duration and Saccade Distance';
    'Fixation and Saccade Statistics by Fixation or Saccade Number';
    '2-D Fixation PDF'
    'Average-Smoothed Fixation Profile';
    'Average-Smoothed Saccade Profile';
    'Persistence Profile';
    'Saccade and Fixation Distance Profiles';
    'Distribution of angles entering a fixation';
    'Distribution of Saccade Amplitudes';
    };

for ff = 1:size(hax,1);
    figure(graphnum+ff)
    set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    for i = 1:4
        switch ff
            case {10,11,12,13}
                if i == 1;
                    pos = [.13 0.57 0.33 0.33];
                elseif i == 2;
                    pos = [.5 0.57 0.33 0.33];
                elseif i == 3;
                    pos = [.13 0.13 0.33 0.33];
                else
                    pos = [.5 0.13 0.33 0.33];
                end
            otherwise
                if i == 1;
                    pos = [.13 0.57 0.33 0.33];
                elseif i == 2;
                    pos = [.57 0.57 0.33 0.33];
                elseif i == 3;
                    pos = [.13 0.13 0.33 0.33];
                else
                    pos = [.57 0.13 0.33 0.33];
                end
        end
        set(hax(ff,i),'position',pos);
    end
    h = findobj(gcf,'Type','axes','Tag','legend');
    if ~isempty(h)
        legpos = get(h,'Position');
        legpos(1) = 0.875;
        legpos(2) = 0.75;
        set(h,'Position',legpos)
    end
    saveas(gcf,[figdir figuretitles{ff}])
    print(gcf,'-r300','-djpeg',[figdir figuretitles{ff}])
end

allviewvariables = {
    'tags: subject names';
    '';
    'allview: all combined data by subject';
    'allview.densitymap: positions of fixations';
    'allview.allfixations: fixation profile by parameters warped twice to median length';
    'allview.allsaccades: saccade profile by parameters warped twice to median length';
    'allview.persistence: persistence/probability of eye movement changing >45 angle double warped';
    'allview.anglebtwfix: angles between fixations by fixation';
    'allview.sacangle_2fix: saccade angles entering a fixation';
    'allview.distanceprofile: velocity of eye movements for saccade + subsequent fixation';
    'allview.distbtwnfix: distance between fixations by fixation number';
    'allview.fixduration: fixation duration by fixation number';
    'allview.sacangle: angle of saccade leaving a fixation by saccade number';
    'allview.sacdist: distance of saccade by saccade number';
    'allview.sacamplitude: saccade amplitude from end to end';
    'allview.timebtwfix: time between fixations by fixation number';
    };

save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\CombinedViewingBehavior.mat'],'tags','allview',...
    'allStatsbyfixation','allviewvariables');
%%
%---[7] Calculate Salience at Returned Fixations---%
% Combines Data by Monkey for Salience at Returned Fixations---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
imageX = 800; imageY = 600;
distancethreshold = [0  24 48 72 96  120  144 168 200 400 800;...%in pixels 24 pixels/dva
                     24 48 72 96 120 144  168 200 400 800 10000];%in pixels 24 pixels/dva
tags = {'MP','TT','JN','IW'};
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'fixation'))
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    eyedatafiles(ii) = i;
                end
            end
        end
    end
    for eyefile = eyedatafiles;
        SalienceIOR(matfiles.mat{eyefile},distancethreshold,imageX,imageY)
    end
end

totalfixations = zeros(1,size(distancethreshold,2));
allSalIOR = cell(1,size(distancethreshold,2));
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'SalienceIOR'))
            load(matfiles.mat{i});
            for ii = 1:length(allSalIOR);
                allSalIOR{ii} = [allSalIOR{ii}; returnfixsal{ii}];
                totalfixations(ii) = totalfixations(ii)+size(returnfixsal{ii},1);
            end
        end
    end
end

%all SalIOR is organized as follows
%fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist fix1dur fix2dur fixnum1 fixnum2
load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'CombinedSalienceStatistics.mat'],'allstatistics','alldata');
chance = allstatistics.confidenceintervals(1,2,1);
allsal = alldata(:,1,1,:);
allsal(isnan(allsal)) = [];
clear alldata allstatistics 

return_sal_pvalues = NaN(size(distancethreshold,2),3);
return_sal_means = NaN(size(distancethreshold,2),2);
return_sal_stds = NaN(size(distancethreshold,2),2);
numreturns = NaN(1,size(distancethreshold,2));
labels = {};
for i = 1:size(distancethreshold,2);
    [~,p]=ttest2(allSalIOR{i}(:,4),allSalIOR{i}(:,8));
    return_sal_pvalues(i,1) = p; %prior vs return
    [~,p]=ttest2(allsal,allSalIOR{i}(:,4));
    return_sal_pvalues(i,2) = p; %prior vs all
    [~,p]=ttest2(allsal,allSalIOR{i}(:,8));
    return_sal_pvalues(i,3) = p; %return vs all
    
    return_sal_means(i,1) = mean(allSalIOR{i}(:,4));
    return_sal_means(i,2) = mean(allSalIOR{i}(:,8));
    return_sal_stds(i,1) = std(allSalIOR{i}(:,4));
    return_sal_stds(i,2) = std(allSalIOR{i}(:,8));
    numreturns(i) = length(allSalIOR{i}(:,8));
    labels{i} = [num2str(distancethreshold(1,i)) '-' num2str(distancethreshold(2,i))];
end

figure
hold on
plot([0 size(distancethreshold,2)+1],[mean(allsal) mean(allsal)],'-k','linewidth',2);
plot([0 size(distancethreshold,2)+1],[mean(chance) mean(chance)],'--k','linewidth',2);
errorbar(return_sal_means(:,2),return_sal_stds(:,2)./sqrt(numreturns'))
errorbar(return_sal_means(:,1),return_sal_stds(:,1)./sqrt(numreturns'),'g')
ylim([0.15 0.45])
for i = 1:size(return_sal_pvalues,1)
    if return_sal_pvalues(i,1) < 0.05 %if the return is
        plot(i,0.43,'r*')
    end
end
xlabel('Distance Tresholds')
ylabel('Salience')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels)
hold off
legend('Mean Salience @ Fixations','Chance Salience','Salience @ Return Fixation',...
    'Salience @ Prior Fixation','p < 0.05','location','northeastoutside')
hold off

load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'CombinedViewingBehavior.mat'],'allview');
fixationduration = [];
for i = 1:length(tags);
    fixationduration = [fixationduration; 5*allview{i}.fixduration];
end
clear allview
fixationduration(isnan(fixationduration)) = [];

return_dur_pvalues = NaN(size(distancethreshold,2),3);
return_dur_means = NaN(size(distancethreshold,2),2);
return_dur_stds = NaN(size(distancethreshold,2),2);
numreturns = NaN(1,size(distancethreshold,2));
labels = {};
for i = 1:size(distancethreshold,2);
    [~,p]=ttest2(allSalIOR{i}(:,10),5*allSalIOR{i}(:,11));
    return_dur_pvalues(i,1) = p; %prior vs return
    [~,p]=ttest2(fixationduration,5*allSalIOR{i}(:,10));
    return_dur_pvalues(i,2) = p; %all vs prior
    [~,p]=ttest2(fixationduration,5*allSalIOR{i}(:,11));
    return_dur_pvalues(i,3) = p; %all vs return
    
    return_dur_means(i,1) = mean(5*allSalIOR{i}(:,10));
    return_dur_means(i,2) = mean(5*allSalIOR{i}(:,11));
    return_dur_stds(i,1) = std(5*allSalIOR{i}(:,10));
    return_dur_stds(i,2) = std(5*allSalIOR{i}(:,11));
    numreturns(i) = length(5*allSalIOR{i}(:,11));
    labels{i} = [num2str(distancethreshold(1,i)) '-' num2str(distancethreshold(2,i))];
end

figure
hold on
plot([0 size(distancethreshold,2)+1],[mean(fixationduration) mean(fixationduration)],'-k','linewidth',1);
errorbar(return_dur_means(:,2),return_dur_stds(:,2)./sqrt(numreturns'))
errorbar(return_dur_means(:,1),return_dur_stds(:,1)./sqrt(numreturns'),'g')

for i = 1:size(return_dur_pvalues,1)
    if return_dur_pvalues(i,1) < 0.05 %if the return is
        plot(i,290,'r*')
    end
end
xlabel('Distance Tresholds')
ylabel('Fixation Duration (ms)')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels)
hold off
legend('Mean Fixation Duration','Duration of Return Fixation',...
    'Duration of Prior Fixation','p < 0.05','location','northeastoutside')
hold off

dur = [];
sal = [];
for i = 1:2
    for ii = 1:size(allSalIOR{i})
        dur = [dur ; 5*allSalIOR{i}(ii,11)]; %prior 10, return 11
        sal = [sal; allSalIOR{i}(ii,8)]; %prior 4, return 8
    end
end
figure
plot(dur,sal,'.')
hold off
xlabel('Duration (ms)')
ylabel('Normalized Salience')
[r,p]= corrcoef(dur,sal);
xlim([0 660])

dist = [];
sal = [];
for i = 1:2
    for ii = 1:size(allSalIOR{i})
        dist = [dist ; allSalIOR{i}(ii,9)];
        sal = [sal; allSalIOR{i}(ii,8)]; %4 for prior, and 8 for return
    end
end
P = polyfit(dist,sal,1);
figure
hold on
plot(dist,sal,'.')
plot(1:50,[1:50]*P(1)+P(2),'r')
hold off
xlabel('Distance (pixels)')
ylabel('Normalized Salience')
[r,p]= corrcoef(dist,sal);

mean_time_btwn_return = [];
mean_numfixations_btwn_return = [];
std_time_btwn_return = [];
std_numfixations_btwn_return = [];
for i = 1:size(distancethreshold,2);
    mean_time_btwn_return(i) = 5*mean(allSalIOR{i}(:,7)-allSalIOR{i}(:,3)-1);
    mean_numfixations_btwn_return(i) = mean(allSalIOR{i}(:,13)-allSalIOR{i}(:,12)-1);
    std_time_btwn_return(i) = 5*std(allSalIOR{i}(:,7)-allSalIOR{i}(:,13)-1);
    std_numfixations_btwn_return(i) = std(allSalIOR{i}(:,13)-allSalIOR{i}(:,12)-1);
end

figure
errorbar(mean_time_btwn_return,std_time_btwn_return./sqrt(numreturns));
xlabel('Distance Tresholds')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels)
ylabel('Duration between prior and return fixation')

figure
errorbar(mean_numfixations_btwn_return,std_numfixations_btwn_return./sqrt(numreturns))
xlabel('Distance Tresholds')
set(gca,'XTick',1:length(distancethreshold))
set(gca,'XTickLabel',labels)
ylabel('Number of fixations in-between prior and return fixation')

returnstats.sal_at_allfixation = allsal;
returnstats.chancesalience = chance;
returnstats.sal_means = return_sal_means;
returnstats.sal_stds = return_sal_stds;
returnstats.sal_pvalues = return_sal_pvalues;
returnstats.dist_sal_correlation = [r(2),p(2)];

returnstats.dur_at_allfixation = fixationduration;
returnstats.dur_means = return_dur_means;
returnstats.dur_stds = return_dur_stds;
returnstats.dur_pvalues = return_dur_pvalues;

returnstats.mean_time_btwn_return = mean_time_btwn_return;
returnstats.std_time_btwn_return = std_time_btwn_return;
returnstats.mean_numfixations_btwn_return = mean_numfixations_btwn_return;
returnstats.std_numfixations_btwn_return = std_numfixations_btwn_return;

SalIORvariablenames = {
    'allSalIOR: combined data from all Salience IOR runs organized by distance';
    'threshold in the following way-[fix1x fix1y fix1t fix1sal fix2x fix2y fix2t...';
    'fix2sal fixdist fix1dur fix2dur fixnum1 fixnum2]';
    '';
    'distancethreshold: distance thresholds to define pairs of fixations';
    '';
    'returnstats....'
    'sal_at_allfixation: salience at fixations 1-35, the median # of fixations';
    'chancesalience: 95% confidence interval of chance salience';
    'sal_means: mean salience by distancthreshold arranged with prior fixations';
    'in the first column and return fixations in second column';
    'sal_stds: std of salience by distancthreshold arranged with prior fixations';
    'in the first column and return fixations in second column';
    'sal_pvalues: p-values comparing salience by distancthreshold arranged with';
    'a kstest comparing salience @ prior fixations to return fixations, salience';
    'at prior fixations to all fixations, and salience at return fixations to all fixations';
    '';
    'dist_sal_correlation: correlation coefficient r, then pvalue of coefficient';
    'dur_:arranged and calculated in the same manner as for salience except for';
    'the relationship between fixation duration and salinece';
    '';
    'Mean_time and number of fixation in between returns and prior fixations'
    };

save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'CombinedSalienceIOR.mat'],'returnstats','totalfixations','allSalIOR',...
    'SalIORvariablenames','distancethreshold')
%%
%---[8] Run BCRW ---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
Combinedbehaviorfile = ['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets'...
    '\CombinedViewingBehavior.mat'];
load(Combinedbehaviorfile,'allview')
imageX = 800; imageY = 600;
plotoptions.runs = 'none'; %all/none
plotoptions.probdens = 'none';
plotoptions.type = 'sal'; %sal/image
% IOR_taus = [1 1/3 1/7 1/12 1/17 1/25 1/35 1/50 0];
IOR_taus = 1/17;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    
    matfiles = what;
    saliencemapfiles = [NaN;NaN];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'saliencemap.mat');
        if ~isempty(str)
            dash = strfind(matfiles.mat{i},'-');
            if ~isempty(str2num(matfiles.mat{i}(dash(1)-1)))
                saliencemapfiles = [saliencemapfiles [i;str2num(matfiles.mat{i}(1:dash(1)-1))]];
            end
        end
    end
    saliencemapfiles(:,1) = [];
    [~,si] = sort(saliencemapfiles(2,:));
    saliencemapfiles = saliencemapfiles(1,si);
    
    for i = 1:length(saliencemapfiles)
        for t = 1:length(tags)
            for it = 1:length(IOR_taus)
                disp(['Running ' tags{t} ' on image #' num2str(i) ' from ' image_sets{SET} ' IOR_tau ' num2str(IOR_taus(it))])
                run_BCRWCF_saveXY(allview{t},matfiles.mat{saliencemapfiles(i)},tags{t},imageX,imageY,plotoptions,IOR_taus(it)) 
                run_CRWCF_saveXY(allview{t},matfiles.mat{saliencemapfiles(i)},tags{t},imageX,imageY,plotoptions,IOR_taus(it)) % same code but no saliene bias
            end
        end
    end
end
%% [9] Calculate goodness of fit of BCRW for fixation location using KL-Divergence
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
binsize = 25; %24 pixels per dva but using 25 cuz divides nicely into 600 & 800
f = fspecial('gaussian',[256,256],24);
KLshuff = NaN(36*length(image_sets),4);
KLnorm = NaN(36*length(image_sets),4);
for imset = 1:length(image_sets);
    disp(['Image set-' num2str(image_sets{imset})])
    dirName = [scm_image_dir image_sets{imset}];
    cd(dirName)
    matfiles = what;
    
    eyedatafiles = zeros(1,length(tags));
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'fixation'))
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    eyedatafiles(ii) = i;
                end
            end
        end
    end
    
    for i = 1:36;
        allfixations = zeros(imageY,imageX);
        allBCRW = zeros(imageY,imageX);
        index = (imset-1)*36+i;
        
        load([num2str(i) '-saliencemap.mat'])
        allsalience = fullmap;
        
        img = double(rgb2gray(imread([num2str(i) '.bmp']))); %want values from 1-256
        imageintensities = 256-img;
        %         imageintensities = img + 1; % if want to use regular image intensities
        
        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;
                
                load(['BCRW IOR TAU 17\' tags{t} '-' num2str(i) '-BCRW.mat'],'fixations')
                allBCRW = allBCRW+fixations;
                
                load(matfiles.mat{eyedatafiles(t)})
                fixations = fixationstats{i*2-1}.fixations;
                if ~isempty(fixations)
                    if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                            fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                        fixations(:,1) = [];
                    end
                    for iii = 1:size(fixations,2)
                        xxyy = round(fixations(:,iii));
                        xxyy(2) = imageY-xxyy(2);
                        xxyy(xxyy < 1) = 1;
                        xxyy(2,(xxyy(2) > imageY)) = imageY;
                        xxyy(1,(xxyy(1) > imageX)) = imageX;
                        allfixations(xxyy(2),xxyy(1)) = allfixations(xxyy(2),xxyy(1)) + 1;
                    end
                end
            end
        end
        allfixations = imfilter(allfixations,f);
        allBCRW = imfilter(allBCRW,f);
        
        binfixations = bin2(allfixations,binsize,binsize);
        binsalience = bin2(allsalience,binsize,binsize);
        binBCRW = bin2(allBCRW,binsize,binsize);
        binimg = bin2(imageintensities,binsize,binsize);
        
        binfixations(binfixations == 0) = eps;
        binsalience(binsalience == 0) = eps;
        binBCRW(binBCRW == 0) = eps;
        binimg(binimg == 0) = eps;
        
        binfixations = binfixations/sum(sum(binfixations));
        binsalience = binsalience/sum(sum(binsalience));
        binBCRW = binBCRW/sum(sum(binBCRW));
        binimg = binimg/sum(sum(binimg));
        
        KLnorm(index,1) = sum(sum(log2(binfixations./binsalience).*binfixations))...
            +sum(sum(log2(binsalience./binfixations).*binsalience));
        KLnorm(index,2) = sum(sum(log2(binfixations./binBCRW).*binfixations))...
            +sum(sum(log2(binBCRW./binfixations).*binBCRW));
        KLnorm(index,3) = sum(sum(log2(binfixations./binimg).*binfixations))...
            +sum(sum(log2(binimg./binfixations).*binimg));
        KLnorm(index,4) = sum(sum(log2(binimg./binsalience).*binimg))...
            +sum(sum(log2(binsalience./binimg).*binsalience));
        
        sz = size(binfixations);
        binsalience = binsalience(randperm(numel(binsalience)));
        binsalience = reshape(binsalience,sz);
        binfixations = binfixations(randperm(numel(binfixations)));
        binfixations = reshape(binfixations,sz);
        binimg = binimg(randperm(numel(binimg)));
        binimg = reshape(binimg,sz);
        binBCRW = binBCRW(randperm(numel(binBCRW)));
        binBCRW = reshape(binBCRW,sz);
        
        KLshuff(index,1) = sum(sum(log2(binfixations./binsalience).*binfixations))...
            +sum(sum(log2(binsalience./binfixations).*binsalience));
        KLshuff(index,2) = sum(sum(log2(binfixations./binBCRW).*binfixations))...
            +sum(sum(log2(binBCRW./binfixations).*binBCRW));
        KLshuff(index,3) = sum(sum(log2(binfixations./binimg).*binfixations))...
            +sum(sum(log2(binimg./binfixations).*binimg));
        KLshuff(index,4) = sum(sum(log2(binimg./binsalience).*binimg))...
            +sum(sum(log2(binsalience./binimg).*binsalience));
    end
end

allKL = KLnorm;
allKLshuff = KLshuff;

pvalues = NaN(1,7);
for i = 1:size(allKL,2);
    [~,p] = kstest2(allKL(:,i),allKLshuff(:,i));
    pvalues(i) = p;
end
[~,p] = kstest2(allKL(:,1),allKL(:,2));
pvalues(5) = p;
[~,p] = kstest2(allKL(:,2),allKL(:,3));
pvalues(6) = p;
[~,p] = kstest2(allKL(:,1),allKL(:,3));
pvalues(7) = p;

pos = 2:3:11;
figure
hold on
b(1) = bar(pos-0.5,mean(allKL),'b');
set(b(1),'BarWidth',0.3)
b(2) = bar(pos+0.5,mean(allKLshuff),'r');
set(b(2),'BarWidth',0.3)
errorbar(pos-0.5,mean(allKL),std(allKL)/sqrt(size(allKL,1)),'xk','linewidth',2);
errorbar(pos+0.5,mean(allKLshuff),std(allKLshuff)/sqrt(size(allKL,1)),'xk','linewidth',2);

for i = 1:length(pvalues);
    if i <= size(allKL,2)
        pos2 = pos(i);
        h = 0.25+max(mean(allKL(:,i))+std(allKL(:,i))/sqrt(size(allKL(:,i),1)),...
            mean(allKLshuff(:,i))+std(allKLshuff(:,i))/sqrt(size(allKLshuff(:,i),1)));
    else
        if i == 5;
            pos2 = 3;
            h = 1+mean(allKLshuff(:,1))+std(allKLshuff(:,1))/sqrt(size(allKLshuff(:,1),1));
            plot([2 2],[h-0.1 h-0.05],'k')
            plot([5 5],[h-0.1 h-0.05],'k')
            plot([2 5],[h-0.05 h-0.05],'k')
        elseif i == 6;
            pos2 = 6.5;
            h = 1+max(mean(allKLshuff(:,3))+std(allKLshuff(:,3))/sqrt(size(allKLshuff(:,3),1)),...
                mean(allKLshuff(:,2))+std(allKLshuff(:,2))/sqrt(size(allKLshuff(:,2),1)));
            plot([5 5],[h-0.1 h-0.05],'k')
            plot([8 8],[h-0.1 h-0.05],'k')
            plot([5 8],[h-0.05 h-0.05],'k')
        elseif i == 7;
            pos2 = 5;
            h = 1.5+max(mean(allKLshuff(:,3))+std(allKLshuff(:,3))/sqrt(size(allKLshuff(:,3),1)),...
                mean(allKLshuff(:,2))+std(allKLshuff(:,2))/sqrt(size(allKLshuff(:,2),1)));
            plot([2 2],[h-0.1 h-0.05],'k')
            plot([8 8],[h-0.1 h-0.05],'k')
            plot([2 8],[h-0.05 h-0.05],'k')
        end
    end
    if pvalues(i) < 0.001
        text(pos2,h,['p = ' num2str(pvalues(i),'%1.1e\n')],'HorizontalAlignment','center')
    else
        text(pos2,h,['p = ' num2str(pvalues(i),'%.3f')],'HorizontalAlignment','center')
    end
end
hold off
legend(b,{'Normal','Shuffled'},'location','NorthEastOutside')
xlim([0 13])
ylim([0 13])
ylabel('KL Divergence: D_{KL} (Distance in Bits)')
set(gca,'XTick',pos)
set(gca,'XTickLabel',{'Fixation vs Salience','Fixation vs BCRW',...
    'Fixation vs Image Intensity','Salainece  vs Image Intensity'})

save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\Combined-KL-DiveregenceTest-CorrectedImgI.mat'],'pvalues','allKL',...
    'allKLshuff')
%%
%---[10] Combine BCRW and Calculate Goodness of Fit for Fixation Location-AUC ROC---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
f = fspecial('gaussian',[256,256],24);
ROC = cell(1,length(image_sets));
all_I_at_fixations = [];
all_BCRW_at_fixations = [];
all_sal_at_fixations = [];
figure
hold on
for imset = 1:length(image_sets);
    ROC{imset} = NaN(3,36);
    disp(['Image set-' num2str(image_sets{imset})])
    dirName = [scm_image_dir image_sets{imset}];
    cd(dirName)
    matfiles = what;
    
    eyedatafiles = zeros(1,length(tags));
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'fixation'))
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    eyedatafiles(ii) = i;
                end
            end
        end
    end

    for i = 1:36;
        salience_at_fixations = [];
        salience_at_random = [];
        BCRW_at_fixations = [];
        BCRW_at_random = [];
        I_at_fixations = [];
        I_at_random = [];
        allBCRW = zeros(imageY,imageX);
        index = (imset-1)*36+i;
        
        load([num2str(i) '-saliencemap.mat'])
        allsalience = fullmap;
        
        imageintensities = 255-double(rgb2gray(imread([num2str(i) '.bmp'])));
        imageintensities = imageintensities/max(max(imageintensities));
        
        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;
                load(['BCRW IOR TAU 17\' tags{t} '-' num2str(i) '-BCRW.mat'],'fixations')
                allBCRW = allBCRW+fixations;
            end
        end
        allBCRW = imfilter(allBCRW,f);
        allBCRW = allBCRW - min(min(allBCRW)); %normalize 0 to 1
        allBCRW = allBCRW/max(max(allBCRW));
        
        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;
                
                load(matfiles.mat{eyedatafiles(t)})
                fixations = fixationstats{i*2-1}.fixations;
                if ~isempty(fixations)
                    if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                            fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                        fixations(:,1) = [];
                    end
                    fixsal  = NaN(1,size(fixations,2));
                    sfixsal = NaN(1,size(fixations,2));
                    fixBCRW = NaN(1,size(fixations,2));
                    sfixBCRW = NaN(1,size(fixations,2));
                    fixI = NaN(1,size(fixations,2));
                    sfixI = NaN(1,size(fixations,2));
                    for iii = 1:size(fixations,2)
                        xxyy = round(fixations(:,iii));
                        xxyy(2) = imageY-xxyy(2);
                        xxyy(xxyy < 1) = 1;
                        xxyy(2,(xxyy(2) > imageY)) = imageY;
                        xxyy(1,(xxyy(1) > imageX)) = imageX;
                        fixsal(iii) = allsalience(xxyy(2),xxyy(1));
                        fixI(iii) = imageintensities(xxyy(2),xxyy(1));
                        fixBCRW(iii) = allBCRW(xxyy(2),xxyy(1));
                        ry = randi(600); ry(ry < 1) = 1;
                        rx = randi(800); rx(rx < 1) = 1;
                        sfixsal(iii) = allsalience(ry,rx);
                        sfixI(iii) = imageintensities(ry,rx);
                        sfixBCRW(iii) = allBCRW(ry,rx);
                    end
                    salience_at_fixations = [salience_at_fixations fixsal];
                    salience_at_random = [salience_at_random sfixsal];
                    BCRW_at_fixations = [BCRW_at_fixations fixBCRW];
                    BCRW_at_random = [BCRW_at_random sfixBCRW];
                    I_at_fixations = [I_at_fixations fixI];
                    I_at_random = [I_at_random sfixI];
                end
            end
        end
        all_I_at_fixations = [all_I_at_fixations I_at_fixations];
        all_BCRW_at_fixations = [all_BCRW_at_fixations BCRW_at_fixations];
        all_sal_at_fixations = [all_sal_at_fixations salience_at_fixations];
        len = length(salience_at_fixations);
        thresh = 0:0.01:1;
        TP = NaN(3,length(thresh)); %True positive
        FA = NaN(3,length(thresh)); %False alarm
        for ii = 1:length(thresh)
            TP(1,ii) = sum(salience_at_fixations > thresh(ii))/len;
            FA(1,ii) = sum(salience_at_random > thresh(ii))/len;
            TP(2,ii) = sum(BCRW_at_fixations > thresh(ii))/len;
            FA(2,ii) = sum(BCRW_at_random > thresh(ii))/len;
            TP(3,ii) = sum(I_at_fixations > thresh(ii))/len;
            FA(3,ii) = sum(I_at_random > thresh(ii))/len;
        end
        ROC{imset}(1,i) = -trapz(FA(1,:),TP(1,:));
        ROC{imset}(2,i) = -trapz(FA(2,:),TP(2,:));
        ROC{imset}(3,i) = -trapz(FA(3,:),TP(3,:));
        plot(FA(1,:),TP(1,:),'b')
        plot(FA(2,:),TP(2,:),'g')
        plot(FA(3,:),TP(3,:),'r')
    end
end
plot(thresh,thresh,'k--','linewidth',3)
hold off
legend('Salience','BCRW','Image Intensity','location','Northeastoutside')
xlabel('False Alarm Rate (FP)')
ylabel('True Positive (TP)')

combinedROC = [];
for imset = 1:length(image_sets);
    combinedROC = [combinedROC ROC{imset}];
end

figure
hold on
bar([1 2 3],mean(combinedROC,2));
errorbar(mean(combinedROC,2),std(combinedROC')/sqrt(size(combinedROC,2)),...
    '+k','linewidth',3);
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel',{'Fixation vs Salience','Fixation vs BCRW',...
    'Fixation vs Image Intensity'})
ylabel('AUC of ROC curve (a.u.)')

zpvalues = NaN(1,size(combinedROC,1));
tpvalues = NaN(1,size(combinedROC,1));
for i = 1:size(combinedROC,1)
    [~,p] = ztest(combinedROC(i,:),0.5,std(combinedROC(i,:)));
    zpvalues(i) = p;
    text(i,mean(combinedROC(i,:))+0.05,['p = ' num2str(p,'%1.1e\n')],'HorizontalAlignment','center')
    if i == 1
        [~,p] = ttest2(combinedROC(1,:),combinedROC(2,:));
        tpvalues(i) = p;
        plot([1.1 1.1],[.77 .8],'k')
        plot([1.9 1.9],[.77 .8],'k')
        plot([1.1 1.9],[.8 .8],'k')
        text(1.45,.8,['p = ' num2str(p,'%1.1e\n')],'HorizontalAlignment','center')
    elseif i == 2
        [~,p] = ttest2(combinedROC(2,:),combinedROC(3,:));
        tpvalues(i) = p;
        plot([2.1 2.1],[.77 .8],'k')
        plot([2.9 2.9],[.77 .8],'k')
        plot([2.1 2.9],[.8 .8],'k')
        text(2.45,.8,['p = ' num2str(p,'%1.1e\n')],'HorizontalAlignment','center')
    elseif i == 3
        [~,p] = ttest2(combinedROC(1,:),combinedROC(3,:));
        tpvalues(i) = p;
        plot([1 1],[.83 .86],'k')
        plot([3 3],[.83 .86],'k')
        plot([1 3],[.86 .86],'k')
        text(2,.86,['p = ' num2str(p,'%1.1e\n')],'HorizontalAlignment','center')
    end
end
hold off
ylim([0 0.75])

save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\Combined-ROC-Corrected-ImageI.mat'],'ROC','zpvalues','tpvalues',...
    'all_I_at_fixations','all_BCRW_at_fixations','all_sal_at_fixations')
%%
%---[11] Combine BCRW and Calculate Goodness of Fit-AUC ROC Fixation Order---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
binsize = 25; %code is optimized for binsize of 25 other values may produce fatal errors
f = fspecial('gaussian',[256,256],24);
parkhurst = {'on','rand','IOR','WTA'}; %distamce bias for fixation order in Parkhurst, Law, Niebur (2002)
for pk = 1:length(parkhurst);
    
    %save set avearges
    ROC = cell(1,length(image_sets));
    all_salience_at_fixations = cell(1,length(image_sets));
    all_salience_at_random = cell(1,length(image_sets));
    all_BCRW_at_fixations = cell(1,length(image_sets));
    all_BCRW_at_random = cell(1,length(image_sets));
    all_I_at_fixations = cell(1,length(image_sets));
    all_I_at_random = cell(1,length(image_sets));
    medianlen = NaN(1,length(image_sets));

    for imset = 1:length(image_sets);
        disp(['Image set-' num2str(image_sets{imset})])
        dirName = [scm_image_dir image_sets{imset}];
        cd(dirName)
        matfiles = what;
        

        salience_at_fixations = NaN(36*length(tags),100);
        salience_at_random = NaN(36*length(tags),100);
        BCRW_at_fixations = NaN(36*length(tags),100);
        BCRW_at_random = NaN(36*length(tags),100);
        I_at_fixations = NaN(36*length(tags),100);
        I_at_random = NaN(36*length(tags),100);
        numfixations = NaN(36*length(tags),1);

        
        eyedatafiles = zeros(1,length(tags));
        for i = 1:length(matfiles.mat);
            if ~isempty(strfind(matfiles.mat{i},'fixation.mat'))
                for ii = 1:length(tags);
                    if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                        eyedatafiles(ii) = i;
                    end
                end
            end
        end
        
        for index = 1:36;
            disp(['Running on image #' num2str(index) ' from ' image_sets{imset}])
            load([num2str(index) '-saliencemap.mat']);
            binsal = bin2(fullmap,binsize,binsize);
            binsal(binsal == 0) = min(min(fullmap));
            salorder = NaN(size(binsal));
            tempsal = binsal;
            [rr,cc] = meshgrid(1:imageX/binsize,1:imageY/binsize);
            fixnum = 1;
            i = 12.5;
            j = 16.5;
            while any(any(binsal > 0));
                if strcmpi(parkhurst{pk},'on')
                    tempsal = tempsal.*exp(-((rr-i).^2+(cc-j).^2)/(2*5^2));
                elseif strcmpi(parkhurst{pk},'rand')
                    tempsal = tempsal.*exp(-((rr-i).^2+(cc-j).^2)/(2*randi([1 10])^2));
                elseif strcmpi(parkhurst{pk},'IOR')
                    C = sqrt((rr-j).^2+(cc-i).^2)<=2;
                    tempsal(C == 1) = 0;
                end
                [i,j] = find(tempsal == max(max(tempsal)));
                if max(max(tempsal)) == 0;
                    binsal = zeros(size(binsal));
                    salorder(isnan(salorder)) = fixnum;
                elseif isempty(i)
                    binsal = zeros(size(binsal));
                    salorder(isnan(salorder)) = fixnum;
                    tempsal = binsal;
                else
                    if length(i) > 1;
                        ind = sub2ind(size(tempsal),i,j);
                        [~,ind2] = sort(binsal(ind),'descend');
                        ind = ind(ind2(1));
                        [i,j] = ind2sub(size(tempsal),ind);
                    end
                    salorder(i,j) = fixnum;
                    binsal(i,j) = 0;
                    tempsal = binsal;
                    tempsal(i,j) = 0;
                    fixnum = fixnum + 1;
                end
            end
            
            img = double(rgb2gray(imread([num2str(index) '.bmp'])));
            img = 256-img;
            binI = bin2(img,binsize,binsize);
            Iorder = NaN(size(binI));
            tempI = binI;
            [rr,cc] = meshgrid(1:imageX/binsize,1:imageY/binsize);
            fixnum = 1;
            i = 12.5;
            j = 16.5;
            while any(any(binI > 0));
                if strcmpi(parkhurst{pk},'on')
                    tempI = tempI.*exp(-((rr-i).^2+(cc-j).^2)/(2*5^2));
                elseif strcmpi(parkhurst{pk},'rand')
                    tempI = tempI.*exp(-((rr-i).^2+(cc-j).^2)/(2*randi([1 10])^2));
                elseif strcmpi(parkhurst{pk},'IOR')
                    C = sqrt((rr-j).^2+(cc-i).^2)<=2;
                    tempI(C == 1) = 0;
                end
                [i,j] = find(tempI == max(max(tempI)));
                if max(max(tempI)) == 0;
                    binI = zeros(size(binI));
                    Iorder(isnan(Iorder)) = fixnum;
                elseif isempty(i)
                    binI = zeros(size(binI));
                    Iorder(isnan(Iorder)) = fixnum;
                    tempI = binI;
                else
                    if length(i) > 1;
                        ind = sub2ind(size(tempI),i,j);
                        [~,ind2] = sort(binI(ind),'descend');
                        ind = ind(ind2(1));
                        [i,j] = ind2sub(size(tempI),ind);
                    end
                    Iorder(i,j) = fixnum;
                    binI(i,j) = 0;
                    tempI = binI;
                    tempI(i,j) = 0;
                    fixnum = fixnum + 1;
                end
            end
            
            maxfixations = 0;
            for t = 1:length(tags)
                load(['BCRW IOR TAU 35\' tags{t} '-' num2str(index) '-BCRW.mat'],'fixationtimes');
                maxfixations = max(maxfixations,max(sum(fixationtimes(:,:,1) > 0,2)));
            end
            fixorderpdf = zeros(imageY,imageX,maxfixations);
            allBCRW = zeros(imageY,imageX);
            for t = 1:length(tags)
                load(['BCRW IOR TAU 35\' tags{t} '-' num2str(index) '-BCRW.mat']);
                allBCRW = allBCRW + fixations;
                for i = 1:size(fixationtimes,1);
                    tind = find(fixationtimes(i,:,1) > 0);
                    for ii = 1:length(tind)
                        x = fixationtimes(i,tind(ii),1);
                        y = fixationtimes(i,tind(ii),2);
                        fixorderpdf(y,x,ii) =  fixorderpdf(y,x,ii) + 1;
                    end
                end
            end
            allBCRW = imfilter(allBCRW,f);
            binallBCRW = bin2(allBCRW,binsize,binsize);
            binfixorderpdf = zeros(imageY/binsize,imageX/binsize,maxfixations);
            for i = 1:maxfixations;
                binfixorderpdf(:,:,i) = bin2(fixorderpdf(:,:,i),binsize,binsize);
            end
            BCRWorder = NaN(imageY/binsize,imageX/binsize);
            for ii = 1:maxfixations
                [i,j,k] = ind2sub(size(binfixorderpdf),find(binfixorderpdf == max(max(max(binfixorderpdf)))));
                if length(k) > 1
                    mind = find(k == min(k));
                    i = i(mind);
                    j = j(mind);
                    k = k(mind);
                    if length(i) > 1
                        ind = sub2ind(size(binallBCRW),i,j);
                        [~,ind2] = sort(binallBCRW(ind),'descend');
                        ind = ind(ind2(1));
                        [i,j] = ind2sub([imageY/binsize,imageX/binsize],ind);
                        k = k(1);
                    end
                end
                binfixorderpdf(:,:,k) = 0;
                binfixorderpdf(i,j,:) = 0;
                BCRWorder(i,j) = k;
            end
            binallBCRW(~isnan(BCRWorder)) = NaN;
            tempBCRW = binallBCRW;
            [rr,cc] = meshgrid(1:imageX/binsize,1:imageY/binsize);
            fixnum = maxfixations+1;
            i = 12.5;
            j = 16.5;
            while any(any(binallBCRW > 0));
                if strcmpi(parkhurst{pk},'on')
                    tempBCRW = tempBCRW.*exp(-((rr-i).^2+(cc-j).^2)/(2*5^2));
                elseif strcmpi(parkhurst{pk},'rand')
                    tempBCRW = tempBCRW.*exp(-((rr-i).^2+(cc-j).^2)/(2*randi([1 10])^2));
                elseif strcmpi(parkhurst{pk},'IOR')
                    C = sqrt((rr-j).^2+(cc-i).^2)<=2;
                    tempBCRW(C == 1) = NaN;
                end
                [i,j] = find(tempBCRW == max(max(tempBCRW)));
                if isnan(min(min(tempBCRW)));
                    binallBCRW(isnan(BCRWorder)) = 0;
                    BCRWorder(isnan(BCRWorder)) = fixnum;
                elseif isempty(i)
                    binallBCRW(isnan(BCRWorder)) = 0;
                    BCRWorder(isnan(BCRWorder)) = fixnum;
                    tempBCRW = binallBCRW;
                else
                    if length(i) > 1;
                        ind = sub2ind(size(tempBCRW),i,j);
                        [~,ind2] = sort(binallBCRW(ind),'descend');
                        ind = ind(ind2(1));
                        [i,j] = ind2sub(size(tempBCRW),ind);
                    end
                    BCRWorder(i,j) = fixnum;
                    binallBCRW(i,j) = NaN;
                    tempBCRW = binallBCRW;
                    tempBCRW(i,j) = 0;
                    fixnum = fixnum + 1;
                end
            end
            
            for t = 1:length(tags)
                disp(['Running ' tags{t} ' on image #' num2str(index) ' from ' image_sets{imset}])
                if eyedatafiles(t) ~= 0;
                    load(matfiles.mat{eyedatafiles(t)})
                    fixations = fixationstats{index*2-1}.fixations;
                    if ~isempty(fixations)
                        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                            fixations(:,1) = [];
                        end
                        
                        loopindex = 36*(t-1)+index;
                        numfixations(loopindex) = size(fixations,2);
                        for iii = 1:size(fixations,2)
                            xxyy = fixations(:,iii);
                            xxyy(2) = imageY-xxyy(2);
                            xxyy = round((xxyy-binsize/2)/binsize+1);
                            xxyy(xxyy < 1) = 1;
                            xxyy(2,(xxyy(2) > size(binI,1))) = size(binI,1);
                            xxyy(1,(xxyy(1) > size(binI,2))) = size(binI,2);
                            
                            salience_at_fixations(loopindex,iii) = abs(salorder(xxyy(2),xxyy(1))-iii);
                            salience_at_random(loopindex,iii) =  abs(salorder(randi(numel(salorder)))-iii);
                            BCRW_at_fixations(loopindex,iii) = abs(BCRWorder(xxyy(2),xxyy(1))-iii);
                            BCRW_at_random(loopindex,iii) = abs(BCRWorder(randi(numel(salorder)))-iii);
                            I_at_fixations(loopindex,iii) = abs(Iorder(xxyy(2),xxyy(1))-iii);
                            I_at_random(loopindex,iii) = abs(Iorder(randi(numel(salorder)))-iii);
                        end
                    end
                end
            end
        end
        
        nans = find(isnan(nanmean(salience_at_fixations)));
        salience_at_fixations(:,nans) = [];
        salience_at_random(:,nans) = [];
        BCRW_at_fixations(:,nans) = [];
        BCRW_at_random(:,nans) = [];
        I_at_fixations(:,nans) = [];
        I_at_random(:,nans) = [];
        
        thresh = 0:1:numel(BCRWorder);
        TP = NaN(3,length(thresh)); %True positive
        FA = NaN(3,length(thresh)); %False alarm
        for fixnum = 1:size(salience_at_fixations,2);
            for ii = 1:length(thresh)
                len = sum(~isnan(salience_at_fixations(:,fixnum)));
                TP(1,ii) = sum(salience_at_fixations(:,fixnum) < thresh(ii))/len;
                FA(1,ii) = sum(salience_at_random(:,fixnum) < thresh(ii))/len;
                TP(2,ii) = sum(BCRW_at_fixations(:,fixnum) < thresh(ii))/len;
                FA(2,ii) = sum(BCRW_at_random(:,fixnum) < thresh(ii))/len;
                TP(3,ii) = sum(I_at_fixations(:,fixnum) < thresh(ii))/len;
                FA(3,ii) = sum(I_at_random(:,fixnum) < thresh(ii))/len;
            end
            ROC{imset}(1,fixnum) = trapz(FA(1,:),TP(1,:));
            ROC{imset}(2,fixnum) = trapz(FA(2,:),TP(2,:));
            ROC{imset}(3,fixnum) = trapz(FA(3,:),TP(3,:));
        end
        
        all_salience_at_fixations{imset} = nanmean(salience_at_fixations);
        all_salience_at_random{imset} = nanmean(salience_at_random);
        all_BCRW_at_fixations{imset} = nanmean(BCRW_at_fixations);
        all_BCRW_at_random{imset} = nanmean(BCRW_at_random);
        all_I_at_fixations{imset} = nanmean(I_at_fixations);
        all_I_at_random{imset} = nanmean(I_at_random);
        medianlen(imset) = nanmedian(numfixations);
        
    end

    medianlen = round(median(medianlen));
    
    figure
    hold on
    for imset = 1:length(image_sets)
        plot(1:medianlen,ROC{imset}(1,1:medianlen),'b')
        plot(1:medianlen,ROC{imset}(2,1:medianlen),'g')
        plot(1:medianlen,ROC{imset}(3,1:medianlen),'r')
    end
    hold off
    xlabel('Fixation Number')
    ylabel('AUC ROC (a.u.)')
    legend('Salience','BCRW','Image Intensity','location','Northeastoutside')
    xlim([0 medianlen])
    title(['Parkurst: ' parkhurst{pk}])
    
    
    figure
    hold on
    for imset = 1:length(image_sets)
        plot(1:medianlen,all_salience_at_fixations{imset}(1:medianlen),'b')
        plot(1:medianlen,all_BCRW_at_fixations{imset}(1:medianlen),'g')
        plot(1:medianlen,all_I_at_fixations{imset}(1:medianlen),'r')
    end
    hold off
    xlabel('Fixation Number')
    ylabel('Difference in Predicted and Actual fixation number')
    legend('Salience','BCRW','Image Intensity','location','Northeastoutside')
    xlim([0 medianlen])
    title(['Parkurst: ' parkhurst{pk}])
    
    if strcmpi(parkhurst{pk},'on')
        save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
            'SCM Image Sets\Combined-FixationOrder-parkhurst-Corrected-ImageI'],...
            'ROC','all_I_at_fixations','all_BCRW_at_fixations',...
            'all_salience_at_fixations','all_I_at_random','all_BCRW_at_random',...
            'all_salience_at_random','medianlen')
    elseif strcmpi(parkhurst{pk},'rand')
        save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
            'SCM Image Sets\Combined-FixationOrder-randparkhurst-Corrected-ImageI'],...
            'ROC','all_I_at_fixations','all_BCRW_at_fixations',...
            'all_salience_at_fixations','all_I_at_random','all_BCRW_at_random',...
            'all_salience_at_random','medianlen')
    elseif strcmpi(parkhurst{pk},'IOR')
        save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
            'SCM Image Sets\Combined-FixationOrder-IOR-Corrected-ImageI'],...
            'ROC','all_I_at_fixations','all_BCRW_at_fixations',...
            'all_salience_at_fixations','all_I_at_random','all_BCRW_at_random',...
            'all_salience_at_random','medianlen')
    elseif strcmpi(parkhurst{pk},'WTA')
        save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
            'SCM Image Sets\Combined-FixationOrder-WTA-Corrected-ImageI'],...
            'ROC','all_I_at_fixations','all_BCRW_at_fixations',...
            'all_salience_at_fixations','all_I_at_random','all_BCRW_at_random',...
            'all_salience_at_random','medianlen')
    else
        error('Method not recognized')
    end
end
%%
% [11.5] Continuation of 11 just sometimes it's easier to run seperately
%---Import data across methods---%
medianlen = 35;
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};

ROCS = cell(4,3);
values = cell(4,3);
medianlen = median(medianlen); %35
load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\Combined-FixationOrder-parkhurst-Corrected-ImageI']);
for ii = 1:length(image_sets)
    for i = 1:3;
        ROCS{1,i} = [ROCS{1,i}; ROC{ii}(i,1:medianlen)];
    end
    values{1,1} = [values{1,1}; all_salience_at_fixations{ii}(1:medianlen)];
    values{1,2} = [values{1,2}; all_BCRW_at_fixations{ii}(1:medianlen)];
    values{1,3} = [values{1,3}; all_I_at_fixations{ii}(1:medianlen)];
end

load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\Combined-FixationOrder-randparkhurst-Corrected-ImageI']);
for ii = 1:length(image_sets)
    for i = 1:3;
        ROCS{2,i} = [ROCS{2,i}; ROC{ii}(i,1:medianlen)];
    end
    values{2,1} = [values{2,1}; all_salience_at_fixations{ii}(1:medianlen)];
    values{2,2} = [values{2,2}; all_BCRW_at_fixations{ii}(1:medianlen)];
    values{2,3} = [values{2,3}; all_I_at_fixations{ii}(1:medianlen)];
end

load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\Combined-FixationOrder-IOR-Corrected-ImageI']);
for ii = 1:length(image_sets)
    for i = 1:3;
        ROCS{3,i} = [ROCS{3,i}; ROC{ii}(i,1:medianlen)];
    end
    values{3,1} = [values{3,1}; all_salience_at_fixations{ii}(1:medianlen)];
    values{3,2} = [values{3,2}; all_BCRW_at_fixations{ii}(1:medianlen)];
    values{3,3} = [values{3,3}; all_I_at_fixations{ii}(1:medianlen)];
end

load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\Combined-FixationOrder-WTA-Corrected-ImageI']);
for ii = 1:length(image_sets)
    for i = 1:3;
        ROCS{4,i} = [ROCS{4,i}; ROC{ii}(i,1:medianlen)];
    end
    values{4,1} = [values{4,1}; all_salience_at_fixations{ii}(1:medianlen)];
    values{4,2} = [values{4,2}; all_BCRW_at_fixations{ii}(1:medianlen)];
    values{4,3} = [values{4,3}; all_I_at_fixations{ii}(1:medianlen)];
end



%---Do an ANOVA analysis to determine what are sources of variance---%
setnum = zeros(length(image_sets),medianlen);
fixnum = zeros(length(image_sets),medianlen);
method = ones(length(image_sets),medianlen);

for i = 1:length(image_sets)
    for ii = 1:medianlen
        setnum(i,ii) = i;
        fixnum(i,ii) = ii;
    end
end
setnum = [setnum(1:end)'; setnum(1:end)'; setnum(1:end)'; setnum(1:end)'];
fixnum = [fixnum(1:end)'; fixnum(1:end)'; fixnum(1:end)'; fixnum(1:end)'];
method = [method(1:end)'; 2*method(1:end)'; 3*method(1:end)'; 4*method(1:end)'];
group = {setnum,fixnum,method};

%can look at interaction too, add "'model','interaction'",but very slow
method_sal_pvalue = anovan(...
    [values{1,1}(1:end)';values{2,1}(1:end)';values{3,1}(1:end)';values{4,1}(1:end)'],...
    group);
method_BCRW_pvalue = anovan(...
    [values{1,2}(1:end)';values{2,2}(1:end)';values{3,2}(1:end)';values{4,2}(1:end)'],...
    group);
method_I_pvalue = anovan(...
    [values{1,3}(1:end)';values{2,3}(1:end)';values{3,3}(1:end)';values{4,3}(1:end)'],...
    group);

meanROCS = zeros(size(ROCS));
stdROCS = zeros(size(ROCS));
for i = 1:size(ROCS,1)
    for ii = 1:size(ROCS,2);
        meanROCS(i,ii) = mean2(ROCS{i,ii});
        stdROCS(i,ii) = std2(ROCS{i,ii});
    end
end

%---Compare AUROCs across maps---%
%compare best methods only
type_pvalues = NaN(1,3);
best(1) = find(meanROCS(:,1) == max(meanROCS(:,1)));
best(2) = find(meanROCS(:,2) == max(meanROCS(:,2)));
best(3) = find(meanROCS(:,3) == max(meanROCS(:,3)));
[~,p] = kstest2(values{best(1),1}(1:end),values{best(2),2}(1:end));
type_pvalues(1) = p;
[~,p] = kstest2(values{best(2),2}(1:end),values{best(3),3}(1:end));
type_pvalues(2) = p;
[~,p] = kstest2(values{best(1),1}(1:end),values{best(3),3}(1:end));
type_pvalues(3) = p;

ROC_pvalues = NaN(1,3);
[~,p] = kstest2(ROCS{best(1),1}(1:end),ROCS{best(2),2}(1:end));
ROC_pvalues(1) = p;
[~,p] = kstest2(ROCS{best(2),2}(1:end),ROCS{best(3),3}(1:end));
ROC_pvalues(2) = p;
[~,p] = kstest2(ROCS{best(1),1}(1:end),ROCS{best(3),3}(1:end));
ROC_pvalues(3) = p;

clr = 'bgr';
fixnum = size(values{1,1},2);
style = {'--','-.','-',':'};
legendlabels = {};
type = {'Salience','BCRW','Image Intensity'};
method = {'Parkhurst 5 dva','Parkhurst 1-10 dva','2 dva IOR Area','Winner Take All'};
figure
hold on
for i = 1:size(values,1) %method
    for ii = 1:size(values,2) %type
        p(i,ii) = plot(1:fixnum,mean(values{i,ii}),[clr(ii) style{i}]);
        legendlabels{ii,i} = [type{ii} ' ' method{i}];
    end
end
hold off
xlim([0 median(medianlen)])
legend(legendlabels{1:end},'location','northeastoutside')
xlabel('Fixation Number')
ylabel('Absolute Difference in Fixation Order')

clr = 'bgr';
fixnum = size(ROCS{1,1},2);
style = {'--','-.','-',':'};
legendlabels = {};
type = {'Salience','BCRW','Image Intensity'};
method = {'Parkhurst 5 dva','Parkhurst 1-10 dva','2 dva IOR Area','Winner Take All'};
figure
hold on
for i = 1:size(ROCS,1) %method
    for ii = 1:size(ROCS,2) %type
        p(i,ii) = plot(1:fixnum,mean(ROCS{i,ii}),[clr(ii) style{i}]);
        legendlabels{ii,i} = [type{ii} ' ' method{i}];
    end
end
hold off
xlim([0 median(medianlen)])
legend(legendlabels{1:end},'location','northeastoutside')
xlabel('Fixation Number')
ylabel('AUC ROC (a.u.)')
ylim([0.5 0.85])

figure
hold on
h = bar(meanROCS');
set(h,'BarWidth',1); % The bars will now touch each other
set(get(gca,'YLabel'),'String','U')
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel',{'Salience','BCRW','Image Intensity'})
ylabel('Difference between actual and predicted fixation order')
legend('Parkhurst 5 dva','Parkhurst 1-10 dva','2 dva IOR Area','Winner Take All',...
    'Location','NortheastOutside')
hold on;
numgroups = size(meanROCS, 2);
numbars = size(meanROCS, 1);
groupwidth = min(0.8, numbars/(numbars+1.5));
numpoints = sqrt(median(medianlen));
for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, meanROCS(i,:), stdROCS(i,:)/numpoints, 'k', 'linestyle', 'none');
end
ylim([0 .85])

text(1.45,0.77,['p = ' num2str(ROC_pvalues(1))],'HorizontalAlignment','center')
plot([1 1],[0.72 0.75],'k')
plot([1.9 1.9],[0.72  0.75],'k')
plot([1 1.9],[0.75 0.75],'k')

text(2.55,0.77,['p = ' num2str(ROC_pvalues(2))],'HorizontalAlignment','center')
plot([2.1 2.1],[0.72 0.75],'k')
plot([3 3],[0.72  0.75],'k')
plot([2.1 3],[0.75 0.75],'k')

text(2,0.825,['p = ' num2str(ROC_pvalues(3))],'HorizontalAlignment','center')
plot([1 1],[0.77 0.8],'k')
plot([3 3],[0.77  0.8],'k')
plot([1 3],[0.8 0.8],'k')