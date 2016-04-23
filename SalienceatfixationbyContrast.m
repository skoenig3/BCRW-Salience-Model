% Determine which type of contrast contributes the most to fixation
% location. Code is derived from fixationSalience_and_significanceCF_contrast.m

scm_image_dir = 'C:\Users\skoenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
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
        fixationSalience_and_significanceCF_contrast(matfiles.mat{eyefile},imageX,imageY,PLOTOPTIONS)
    end
end
clearvars -except scm_image_dir image_sets

scm_image_dir = 'C:\Users\skoenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
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
        str = strfind(matfiles.mat{i},'FixationStatistics_contrast');
        if ~isempty(str)
            statfiles = [statfiles i];
        end
    end
    for stat = statfiles;
        load(matfiles.mat{stat},'statistics')
        minlen = min(minlen,size(statistics.numbervalues,2));
    end
end

alldata = NaN(36*length(image_sets)*4,4,minlen);
allshuffled = NaN(36*length(image_sets)*4,4,minlen);
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    statfiles = zeros(1,length(tags));
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'FixationStatistics'));
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    statfiles(ii) = i;
                end
            end
        end
    end
    
    for stat = 1:4;
        if statfiles(stat) ~= 0
            i = SET*stat;
            load(matfiles.mat{statfiles(stat)})
            combineddata = shuffunshuffdata{2}{3};
            combinedshuffled = shuffunshuffdata{1}{3};
            alldata(36*(i-1)+1:i*36,:,:) = combineddata(:,:,1:minlen);
            allshuffled(36*(i-1)+1:i*36,:,:) = combinedshuffled(:,:,1:minlen);
        end
    end
end

allmeanvals = NaN(4,minlen);
allstdvals = NaN(4,minlen);
allnumvals = NaN(4,minlen);
for i = 1:size(allmeanvals,1)
    for ii = 1:minlen
        allmeanvals(i,ii) = nanmean(alldata(:,i,ii));
        allstdvals(i,ii) = nanstd(alldata(:,i,ii));
        allnumvals(i,ii) = sum(~isnan(alldata(:,i,ii)));
    end
end

% z-test of means agains random distributions assuming mean is larger
allzp = NaN(size(allmeanvals)); %p-values
allcI = NaN(size(allmeanvals)); %top confidence interval value, lowest is typticall 0/-Inf
for i = 1:size(allmeanvals,1)
    for ii = 1:size(allmeanvals,2)
        shuffledvals = allshuffled(:,i,:);
        shuffledvals(isnan(shuffledvals)) = [];
        [~,p,ci] = ztest(shuffledvals,allmeanvals(i,ii),std(shuffledvals),...
            0.05);
        allzp(i,ii) = p;
        allcI(i,ii) = ci(2);
    end
end

allstatistics.meanvalues = allmeanvals;
allstatistics.stdvalues = allstdvals;
allstatistics.numbervalues = allnumvals;
allstatistics.pvalues = allzp;
allstatistics.confidenceintervals = allcI;

legendlabels = {'Salience','Intensity Contrast','Color Contrast','Orientation Contrast'};
averaginglabels = {' at mean fixation location'};
for i = 1:size(allcI,1)
    figure
    hold on
    p = plot(allcI(i,:),'r--');
    pp = errorbar(allmeanvals(i,:),allstdvals(i,:)./sqrt(allnumvals(i,:)),'r');
    hold off
    title([legendlabels{i} ' at a fixation by fixation number across all'...
        ' all monkeys and image sets'])
    legend([p pp],{['Chance ' legendlabels{i} averaginglabels{1}],...
        [legendlabels{i} averaginglabels{1}]},'Location','NorthEastOutside');
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
    'parameters during a fixatoin. Z-column is organized by fixation number';
    };

save(['C:\Users\skoenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'CombinedSalienceStatistics_contrast.mat'],'alldata','allshuffled',...
    'allstatistics','allstatvariablenames');
clearvars -except scm_image_dir image_sets