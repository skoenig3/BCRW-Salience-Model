function fixationSalience_and_significanceCF_contrast(FIXATIONFILE,imageX,imageY,PLOTOPTIONS)
% created by Seth Koenig 6/21/2012 modified for Detecing saccades with Cluster Fixation (CF)
% on 11/15/2012

% function determines normalized salience, salience contrast values, and image
% intensity values at the location of fixation. The function also calculates
% if fixations occur at these values at a probability greater than expected
% by chance using both a z-test and a t-test. Values are plotted as a
% function of fixation number.

% Inputs:
%   FIXATIONFILE: Fixations extracted from cortex e.g. MP120606_1-fixation.mat
%   imageX & imageY: X and Y dimensions of the set's images,respstivley
%
%   PLOTOPTIONS: 'all' if want a plot containing mean and standard deviations
%   normalized salience, salience contrast values, and image intensity values
%   at the location of fixation in comparison to their shuffled counterparts.
%   Stastics are calculated against entire random distribution not against
%   individual fixation numbers (some shuffled data are randomly higher
%   than others).

%Outputs:
%   A .mat file named [FIXATIONFILE(1:end-14) '-FixationStatistics']
%   containg fixations statistics. See variable statvariablenames for
%   detailed explanation of variables in .mat file.

if nargin < 3
    error(['Not enough inputs: function requires ViewingBehaviorFile,'...
        'imageX, imageY, and PLOTOPTIONS.'])
elseif nargin < 4
    PLOTOPTIONS = 'none';
end

load(FIXATIONFILE);

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

maxfixs = 0;
numtrials = round(length(fixationstats)/2);
shuffunshuff = cell(1,2); %1 for shuffled 2 for unshuffled data
for s = 1:length(shuffunshuff)
    corrs = cell(numtrials,4);
    for cndlop=1:2:length(fixationstats)
        reindexed = (cndlop+1)/2;
        load(matfiles.mat{saliencemapfiles(reindexed)});
        
        fixationtimes = fixationstats{cndlop}.fixationtimes;
        if ~isempty(fixationtimes)
            fixations = fixationstats{cndlop}.fixations;
            xy = fixationstats{cndlop}.XY;
            if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                    fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                fixations(:,1) = [];
                fixationtimes(:,1) = [];
            end
            
            numfixs = size(fixations,2);
            maxfixs = max(maxfixs,numfixs);
            for i = 1:size(corrs,2)
                for ii = 1:size(corrs,3)
                    corrs{reindexed,i,ii} = NaN(1,numfixs);
                end
            end
            
            for i = 1:numfixs
                if s == 1; %shuffled points
                    spot = [ceil(imageX*rand) ceil(imageY*rand)]; %fake x,y data
                else %unshuffled
                    spot = ceil(fixations(:,i));
                    spot(2) = imageY-spot(2);
                    spot(spot < 1) = 1;
                    spot(1,spot(1) > imageX) = imageX;
                    spot(2,spot(2) > imageY) = imageY;
                end
                corrs{reindexed,1}(i) = fullmap(spot(2),spot(1));
                corrs{reindexed,2}(i) = Imap(spot(2),spot(1));
                corrs{reindexed,3}(i) = Cmap(spot(2),spot(1));
                corrs{reindexed,4}(i) = Omap(spot(2),spot(1));
            end
        end
    end
    shuffunshuff{s} = corrs;
end

for s = 1:length(shuffunshuff)
    combineddata = NaN(numtrials,4,maxfixs);
    for i = 1:size(combineddata,1)
        for ii = 1:size(combineddata,2)
            if ~isempty(shuffunshuff{s}{i,ii})
            combineddata(i,ii,1:length(shuffunshuff{s}{i,ii}))= ...
                shuffunshuff{s}{i,ii};
            end
        end
    end
    
    meanvals = NaN(4,maxfixs);
    stdvals = NaN(4,maxfixs);
    numvals = NaN(4,maxfixs);
    for i = 1:size(meanvals,1)
        for ii = 1:maxfixs
            if  sum(~isnan(combineddata(:,i,ii))) > 5;
                meanvals(i,ii) = nanmean(combineddata(:,i,ii));
                stdvals(i,ii) = nanstd(combineddata(:,i,ii));
                numvals(i,ii) = sum(~isnan(combineddata(:,i,ii)));
            end
        end
    end
    [i ii] = ind2sub(size(meanvals),find(isnan(meanvals)));
    ii = unique(ii);
    meanvals(:,ii) = []; stdvals(:,ii) = []; numvals(:,ii) = [];
    
    shuffunshuffdata{s} = {meanvals stdvals combineddata};
end

% z-test of means agains random distributions assuming mean is larger
zp = NaN(size(meanvals)); %p-values
cI = NaN(size(meanvals)); %top confidence interval value, lowest is typticall 0/-Inf
for i = 1:size(meanvals,1)
    for ii = 1:size(meanvals,2)
        shuffleddata =  shuffunshuffdata{1}{3}(:,i,ii,:);
        shuffleddata(isnan(shuffleddata)) = [];
        for iii = 1:size(meanvals,3)
            [~,p,ci] = ztest(shuffleddata,meanvals(i,ii,iii),std(shuffleddata),...
                0.05,'left');
            zp(i,ii,iii) = p;
            cI(i,ii,iii) = ci(2);
        end
    end
end

statistics.meanvalues = meanvals;
statistics.stdvalues = stdvals;
statistics.numbervalues = numvals;
statistics.pvalues = zp;
statistics.confidenceintervals = cI;

if strcmpi(PLOTOPTIONS,'all')
    clrs = ['rb'];
    legendlabels = {'Salience','Salience Contrast','Image Intensity'};
    averaginglabels = {' at mean fixation location',' mean during fixation'};
    for i = 1:size(cI,1)
        figure
        hold on
        for ii = 1:size(cI,2)
            p(ii)= plot(cI(i,ii)*ones(1,size(cI,3)),['--' clrs(ii)]);
            pp(ii) = errorbar(shiftdim(meanvals(i,ii,:)),...
                shiftdim(stdvals(i,ii,:))./sqrt(shiftdim(numvals(i,ii,:))),clrs(ii));
        end
        hold off
        title([FIXATIONFILE(1:end-13) 'Image ' legendlabels{i}...
            ' at a fixation by fixation number'])
        legend([p pp],{['Chance ' legendlabels{i} averaginglabels{1}],...
            ['Chance ' legendlabels{i} averaginglabels{2}],...
            [legendlabels{i} averaginglabels{1}], ...
            [legendlabels{i} averaginglabels{2}]},'Location','NorthEastOutside')
        xlabel('Fixation Number')
        ylabel('Normalized Value')
    end
end

statvariablenames = {
    'shuffunshuffdata: cell arrray containing shuffled {1} and unshuffled {2} values.';
    'The cell arrays contain meanvalues, stdvalues and all combined values for';
    'salience, salience contrast, and image intensity at fixations, respectively.';
    'Each of these cells are arranged by row,column,and z-column in the following';
    'manner:  Rows are arranged as Salience, salience contrast, and intensity;';
    'Columns indicate if data is these parmaters at the average fixation cooridante';
    '(col 1) or the average of the parameters during a fixatoin; Z-column is';
    'organized by fixation number';
    '';
    'statistics: structure with results from a z-test to determine if fixations';
    'occur at salience, salience contrasts, and image intesntisy valeus at rates';...
    'higher than what would be expected by chance. Random distributions from each';
    'parameter is compared to the mean value at that parameter by fixation';
    'statistics contains means, std, number of fixations, p-values,and confidence intervals';
    'These variables are arranged by row,column,z-column. Rows are arranged as';
    'Salience, salience contrast, and intensity. Columns indicate if data is these';
    'paramaters at the average fixation cooridante (col 1) or the average of the';
    'parameters during a fixatoin. Z-column is organized by fixation number';
    };

save([FIXATIONFILE(1:end-13) '-FixationStatistics_contrast'],...
    'shuffunshuffdata','statistics','statvariablenames')
end