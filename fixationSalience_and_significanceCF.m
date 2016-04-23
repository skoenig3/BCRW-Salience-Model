function fixationSalience_and_significanceCF(FIXATIONFILE,imageX,imageY,PLOTOPTIONS)
% created by Seth Koenig 6/21/2012 modified for Detecing saccades with Cluster Fixation (CF)
% on 11/15/2012

% function determines normalized salience, salience contrast values, and image
% intensity values at the location of fixation. The function also calculates
% if fixations occur at these values at a probability greater than expected
% by chance using both a z-test and a t-test. Values are plotted as a
% function of fixation number. Salience Contrast is now unused; It was a 
% just an idea I wanted to try out.

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

load(FIXATIONFILE,'fixationstats');

RF = receptivefield(imageY,imageX);%"receptive field" to calclulate salience contrast

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
    corrs = cell(numtrials,3,2);
    for cndlop=1:2:length(fixationstats)
        reindexed = (cndlop+1)/2;
        [saliencemap exe img] = getMaps(matfiles.mat{saliencemapfiles(reindexed)},RF);
        
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
                    L = diff(fixationtimes(:,i));
                    spot = [ceil(imageX*rand) ceil(imageY*rand)]; %fake x,y data
                    spots = [ceil(imageX*rand(1,L));ceil(imageY*rand(1,L))];
                    spots = sub2ind([imageY imageX],spots(2,:),spots(1,:));
                else %unshuffled
                    spot = ceil(fixations(:,i));
                    spot(2) = imageY-spot(2);
                    spot(spot < 1) = 1;
                    spot(1,spot(1) > imageX) = imageX;
                    spot(2,spot(2) > imageY) = imageY;
                    spots = [xy(1,fixationtimes(1,i):fixationtimes(2,i));...
                        xy(2,fixationtimes(1,i):fixationtimes(2,i))];
                    spots(2,:) = imageY-spots(2,:);
                    spots = round(spots);
                    spots(spots < 1) = 1;
                    spots(1,spots(1,:) > imageX) = imageX;
                    spots(2,spots(2,:) > imageY) = imageY;
                    spots = sub2ind([imageY imageX],spots(2,:),spots(1,:));
                end
                corrs{reindexed,1,1}(i) = saliencemap(spot(2),spot(1));
                corrs{reindexed,1,2}(i) = mean(saliencemap(spots));
                corrs{reindexed,2,1}(i) = exe(spot(2),spot(1));
                corrs{reindexed,2,2}(i) = mean(exe(spots));
                corrs{reindexed,3,1}(i) = img(spot(2),spot(1));
                corrs{reindexed,3,2}(i) = mean(img(spots));
            end
        end
    end
    shuffunshuff{s} = corrs;
end

for s = 1:length(shuffunshuff)
    combineddata = NaN(numtrials,3,2,maxfixs);
    for i = 1:size(combineddata,1)
        for ii = 1:size(combineddata,2)
            for iii = 1:size(combineddata,3)
                if ~isempty(shuffunshuff{s}{i,ii,iii})
                    combineddata(i,ii,iii,1:length(shuffunshuff{s}{i,ii,iii}))= ...
                        shuffunshuff{s}{i,ii,iii};
                end
            end
        end
    end
    
    meanvals = NaN(3,2,maxfixs);
    stdvals = NaN(3,2,maxfixs);
    numvals = NaN(3,2,maxfixs);
    for i = 1:size(meanvals,1)
        for ii = 1:size(meanvals,2)
            for iii = 1:maxfixs
                if  sum(~isnan(combineddata(:,i,ii,iii))) > 5;
                    meanvals(i,ii,iii) = nanmean(combineddata(:,i,ii,iii));
                    stdvals(i,ii,iii) = nanstd(combineddata(:,i,ii,iii));
                    numvals(i,ii,iii) = sum(~isnan(combineddata(:,i,ii,iii)));
                end
            end
        end
    end
    [i ii iii] = ind2sub(size(meanvals),find(isnan(meanvals)));
    iii = unique(iii);
    meanvals(:,:,iii) = []; stdvals(:,:,iii) = []; numvals(:,:,iii) = [];
    
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
    '(col 1) or the average of the parameters during a fixation; Z-column is';
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

save([FIXATIONFILE(1:end-13) '-FixationStatistics'],...
    'shuffunshuffdata','statistics','statvariablenames')

    function [RF] = receptivefield(rr,cc)
        % rr, cc unused. use to be be and were imageY, and imageX respectively
        rfsize = 30; %diameter ~2.5 degrees
        sig1 = 5;
        sig2 = 10; %sigma changes not a huge effect on results, RF bigger change
        [x,y] = meshgrid(-rfsize:rfsize);
        gabor1 = exp(-(x.^2 + y.^2)/(2*sig1^2));
        gabor2 = exp(-(x.^2 + y.^2)/(2*sig2^2));
        gabor1 = gabor1/sum(sum(gabor1));
        gabor2 = gabor2/sum(sum(gabor2));
        RF = gabor1-gabor2;
    end

    function [saliencemap,exe,img] = getMaps(salmapfile,RF)
        dash = strfind(salmapfile,'-');
        imgfile = [salmapfile(1:dash(1)-1) '.bmp'];
        img = imread(imgfile);
        img = double(rgb2gray(img));
        img = img - min(min(img)); 
        img= img/max(max(img)); %normalize img from [0 1]
        load(salmapfile,'fullmap');
        saliencemap = fullmap;
        exe = imfilter(fullmap,RF,'replicate'); %salience contrast map
        exe = exe - min(min(exe)); %normalize & scale from [0 1]
        exe = exe/max(max(exe));
    end
end