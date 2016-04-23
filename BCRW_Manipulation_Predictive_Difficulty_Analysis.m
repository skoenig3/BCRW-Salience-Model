% Code written By Seth Konig to determine if certain maniuplations are
% harder than others to spot. Written July 2015

%%
%---[1] Make Salience maps for manipulated imagse---%
% scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
% image_sets = {'Set006','Set007','Set008','Set009',...
%     'SetE001','SetE002','SetE003','SetE004'};
%
%
% type = ['om']; %object replaced and object moved
%
% for SET =1:length(image_sets);
%     cd([scm_image_dir image_sets{SET} '\']);
%     for img = 1:36
%         for t = 1:2
%             getSalienceMap([num2str(img) type(t) '.bmp'])
%         end
%     end
% end
% emailme('Finished creating saliencemaps')
%%
%---[2] Run simulations for manipulated images---%
% scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
% image_sets = {'Set006','Set007','Set008','Set009',...
%     'SetE001','SetE002','SetE003','SetE004'};
% tags = {'MP','TT','JN','IW'};
% Combinedbehaviorfile = ['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets'...
%     '\CombinedViewingBehavior.mat'];
% load(Combinedbehaviorfile,'allview')
% imageX = 800; imageY = 600;
% plotoptions.runs = 'none'; %all/none
% plotoptions.probdens = 'none';
% plotoptions.type = 'sal'; %sal/image
% IOR_tau = [1/17];
%
% %same name for everyset
%
% saliencemapfiles = {};
% type = ['om']; %object replaced and object moved
% for img = 1:36
%     for t = 1:2
%         saliencemapfiles{img,t} = [num2str(img) type(t) '-saliencemap.mat'];
%     end
% end
% saliencemapfiles= saliencemapfiles(1:end);
%
% for SET = 1:length(image_sets);
%     SETNUM = image_sets{SET};
%     cd([scm_image_dir SETNUM])
%
%
%     for i = 1:length(saliencemapfiles)
%         for t = 1:length(tags)
%             disp(['Running ' tags{t} ' on image #' num2str(i) ' from ' image_sets{SET}])
%             run_BCRWCF(allview{t},saliencemapfiles{i},tags{t},imageX,imageY,plotoptions,IOR_tau)
%         end
%     end
% end
%%
%---[3] Determine the number of fixations to Manipulated ROI---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
load('C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\CriticalRegions.mat')%SCM ROIs

type = ['om']; %object replaced and object moved
imageY = 600;
imageX = 800;

replaced_difficulty = NaN(length(image_sets),36);
moved_difficulty = NaN(length(image_sets),36);

replaced_total = NaN(length(image_sets),36);
moved_total = NaN(length(image_sets),36);
replaced_area = NaN(length(image_sets),36);
moved_area = NaN(length(image_sets),36);
%area should be the same since the same object is being manipulated
for SET =1:length(image_sets);
    cd([scm_image_dir image_sets{SET} '\']);
    
    for index = 1:36
        for s = 1:2
            if isnan(allROIs{SET,s}(index,1))
                disp('image skipped')
                continue
            end
            fixations = cell(1,length(tags));
            for t = 1:length(tags)
                load(['BCRW IOR TAU 17\' tags{t} '-' num2str(index) type(s) '-BCRW.mat'],'fixationtimes');
                for i = 1:size(fixationtimes,1);
                    fixations{t}{i} = NaN(2,60);
                end
                for i = 1:size(fixationtimes,1);
                    tind = find(fixationtimes(i,:,1) > 0);
                    if length(tind) > 60
                        tind = tind(1:60);
                    end
                    for ii = 1:length(tind)
                        x = fixationtimes(i,tind(ii),1);
                        y = fixationtimes(i,tind(ii),2);
                        fixations{t}{i}(:,ii) = [x;y];
                    end
                end
            end
            
            difficulty = NaN(4,length(fixations{1}));
            total = NaN(4,length(fixations{1}));
            for t = 1:length(fixations)
                for i = 1:length(fixations{t})
                    %add 1 dva or 24 pixel buffer around ROI
                    ind = find(fixations{t}{i}(1,:)+24 > allROIs{SET,s}(index,1)...
                        & fixations{t}{i}(1,:)-24 < allROIs{SET,s}(index,2) ....
                        & fixations{t}{i}(2,:)+24 > allROIs{SET,s}(index,3) ...
                        & fixations{t}{i}(2,:)-24 < allROIs{SET,s}(index,4));%BCRW y coordinates already flipped
                    if ~isempty(ind)
                        difficulty(t,i) = min(ind);
                        total(t,i) = length(ind);
                    else
                        difficulty(t,i) = NaN; %if not found cap
                        total(t,i) = 0;
                    end
                end
            end
            img = zeros(600,800);
            for t = 1:length(fixations)
                for i = 1:length(fixations{t})
                    for ii = 1:size(fixations{t}{i},2)
                        if ~isnan(fixations{t}{i}(1,ii))
                            img(fixations{t}{i}(2,ii),fixations{t}{i}(1,ii)) =  img(fixations{t}{i}(2,ii),fixations{t}{i}(1,ii))+1;
                        end
                    end
                end
            end
            if s == 1
                replaced_difficulty(SET,index) = nanmean(difficulty(1:end));
                replaced_total(SET,index) = nanmean(total(1:end));
                replaced_area(SET,index) = (allROIs{SET,s}(index,2)-allROIs{SET,s}(index,1))*...
                    (allROIs{SET,s}(index,4)-allROIs{SET,s}(index,3));
            else
                moved_difficulty(SET,index) = nanmean(difficulty(1:end));
                moved_total(SET,index) = nanmean(total(1:end));
                moved_area(SET,index) = (allROIs{SET,s}(index,2)-allROIs{SET,s}(index,1))*...
                    (allROIs{SET,s}(index,4)-allROIs{SET,s}(index,3));
            end
        end
    end
end

nanmean(nanmean(replaced_difficulty))
nanmean(nanmean(moved_difficulty))
nanmean(nanmean(replaced_total))
nanmean(nanmean(moved_total))
nanmean(nanmean(replaced_area))
nanmean(nanmean(moved_area))

[~,tpd] = ttest2(replaced_difficulty(1:end),moved_difficulty(1:end))
[~,kpd] = kstest2(replaced_difficulty(1:end),moved_difficulty(1:end))
[~,tpt] = ttest2(replaced_total(1:end),moved_total(1:end))
[~,kpt] = kstest2(replaced_total(1:end),moved_total(1:end))
[~,tpa] = ttest2(replaced_area(1:end),moved_area(1:end))
[~,kpa] = kstest2(replaced_area(1:end),moved_area(1:end))

%%
% %esther's sets versus wills sets
% wdr = replaced_difficulty(1:4,:);
% edr = replaced_difficulty(5:end,:);
% wdm= moved_difficulty(1:4,:);
% edm = moved_difficulty(5:end,:);
% 
% [~,p] = ttest2(wdr(1:end),edr(1:end))
% [~,p] = kstest2(wdr(1:end),edr(1:end))
% [~,p] = ttest2(wdm(1:end),edm(1:end))
% [~,p] = kstest2(wdm(1:end),edm(1:end))
% 
% wdr = replaced_total(1:4,:);
% edr = replaced_total(5:end,:);
% wdm= moved_total(1:4,:);
% edm = moved_total(5:end,:);
% 
% [~,p] = ttest2(wdr(1:end),edr(1:end))
% [~,p] = kstest2(wdr(1:end),edr(1:end))
% [~,p] = ttest2(wdm(1:end),edm(1:end))
% [~,p] = kstest2(wdm(1:end),edm(1:end))
%%

scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
load('C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\CriticalRegions.mat')%SCM ROIs

type = ['om']; %object replaced and object moved
imageY = 600;
imageX = 800;

moved_fixations = NaN(length(image_sets),36,length(tags)); %proportion total number of fixations
moved_fix10 = NaN(length(image_sets),36,length(tags)); %proportion of of 1st 10 fixations
moved_time = NaN(length(image_sets),36,length(tags)); %proportion of total time
moved_time3 = NaN(length(image_sets),36,length(tags));%proprotion of 1st 3 seconds
moved_time15 = NaN(length(image_sets),36,length(tags));%proprotion of 1st 3 seconds

replaced_fixations = NaN(length(image_sets),36,length(tags));  %proportion total number of fixations
replaced_fix10 = NaN(length(image_sets),36,length(tags));%proportion of of 1st 10 fixations
replaced_time = NaN(length(image_sets),36,length(tags));%proportion of total time
replaced_time3 = NaN(length(image_sets),36,length(tags));%proprotion of 1st 3 seconds

for SET = 1:length(image_sets)
    cd([scm_image_dir image_sets{SET} '\']);
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation.mat');
        if ~isempty(str)
            if SET >= 5 && strcmpi(matfiles.mat{i}(1:2),tags{2}) %TT has lesion for esthers sets
                continue
            end
            eyedatafiles = [eyedatafiles i];
        end
    end
    for eye = 1:length(eyedatafiles)
        load(matfiles.mat{eyedatafiles(eye)})
        
        for i = 1:length(tags)
            if strcmpi(matfiles.mat{eyedatafiles(eye)}(1:2),tags{i})
               monk = i;
               break
            end
        end
        
        for i = 2:2:length(fixationstats);
            if trialtype(i) == 'f' || trialtype(i) == 'e'
                continue
            end
            if isempty(fixationstats{i})
               continue 
            end
            
            %exclude images that do not meet SCM criteria (based on earlier
            %analysis done by Harkirat)
            if SET==1 %Set006
                if i==4
                    continue;
                end
            end
            if SET==2 %Set007
                if i==4
                    continue;
                end
            end
            if SET==3 %Set008
                if i==38
                    continue;
                end
            end
            if SET==4 %Set009
                if i==4 || i==34
                    continue;
                end
            end
            
            novfixations = fixationstats{i-1}.fixations;
            novtimes =  5*fixationstats{i-1}.fixationtimes;
            novfixations(2,:) = imageY-novfixations(2,:);
            
            repfixations = fixationstats{i}.fixations;
            reptimes =  5*fixationstats{i}.fixationtimes;
            repfixations(2,:) = imageY-repfixations(2,:);
            
            novtimes(novtimes > 10000) = 10000;
            reptimes(reptimes > 10000) = 10000;
            
            if any(isnan(allROIs{SET,1}(i/2,:)))
                continue
            end
            if strcmp(trialtype(i),'r')
                ROI1 = allROIs{SET,1}(i/2,:); %[min_c max_c min_r max_r]
                ROI1(1) = ROI1(1)-24;
                ROI1(2) = ROI1(2)+24;
                ROI1(3) = ROI1(3)-24;
                ROI1(4) = ROI1(4)+24;
                ROI1(ROI1 < 1) = 1;
                ROI1(ROI1 > imageX) = imageX;
                if ROI1(3) > imageY; ROI1(3) = imageY;end
                if ROI1(4) > imageY; ROI1(4) = imageY;end
                ROI2 =ROI1;
                type = 1;
            elseif strcmp(trialtype(i),'m');
                ROI1 = allROIs{SET,3}(i/2,:); % before moved ROI
                ROI1(1) = ROI1(1)-24;
                ROI1(2) = ROI1(2)+24;
                ROI1(3) = ROI1(3)-24;
                ROI1(4) = ROI1(4)+24;
                ROI1(ROI1 < 1) = 1;
                ROI1(ROI1 > imageX) = imageX;
                if ROI1(3) > imageY; ROI1(3) = imageY;end
                if ROI1(4) > imageY; ROI1(4) = imageY;end
                ROI2 = allROIs{SET,2}(i/2,:); %moved ROI
                ROI2(1) = ROI2(1)-24;
                ROI2(2) = ROI2(2)+24;
                ROI2(3) = ROI2(3)-24;
                ROI2(4) = ROI2(4)+24;
                ROI2(ROI2 < 1) = 1;
                ROI2(ROI2 > imageX) = imageX;
                if ROI2(3) > imageY; ROI2(3) = imageY;end
                if ROI2(4) > imageY; ROI2(4) = imageY;end
                type = 2;
            end
            
            ROIfix2 = find(...
                repfixations(1,:) > ROI2(1) & repfixations(1,:) < ROI2(2) & ...
                repfixations(2,:) > ROI2(3) & repfixations(2,:) < ROI2(4));%ROIfix2 contains indices of repfixations corresponding to fixations that occur in ROI of the manipulated image
            
            ROIfix1 = find(...
                novfixations(1,:) > ROI1(1) & novfixations(1,:) < ROI1(2) & ...
                novfixations(2,:) > ROI1(3) & novfixations(2,:) < ROI1(4)); %ROIfix1 contains indices of novfixations corresponding to fixations that occur in ROI of novel image.
            
            if isempty(ROIfix1)
                continue
            end
            
            tempvec=zeros(1,10000);
            for j=1:length(ROIfix2)
                tempvec(reptimes(1,ROIfix2(j)):reptimes(2,ROIfix2(j)))=1;
            end
            
            
            if strcmp(trialtype(i),'m');
                moved_fixations(SET,i/2,monk) = length(ROIfix2)/length(repfixations);
                moved_fix10(SET,i/2,monk) =   sum(ROIfix2 < 10)/10;
                moved_time(SET,i/2,monk) = sum(tempvec)/10000;
                moved_time3(SET,i/2,monk) = sum(tempvec(1:3000))/3000;
                moved_time15(SET,i/2,monk) = sum(tempvec(500:1500))/1000;
            else
                replaced_fixations(SET,i/2,monk) = length(ROIfix2)/length(repfixations);
                replaced_fix10(SET,i/2,monk) =   sum(ROIfix2 < 10)/10;
                replaced_time(SET,i/2,monk) = sum(tempvec)/10000;
                replaced_time3(SET,i/2,monk) = sum(tempvec(1:3000))/3000;
            end
        end
    end
end
%%

moved_fixations = nanmean(moved_fixations,3);
moved_fix10 = nanmean(moved_fix10,3);
moved_time = nanmean(moved_time,3);
moved_time3 = nanmean(moved_time3,3);
moved_time15 = nanmean(moved_time15,3);

replaced_fixations = nanmean(replaced_fixations,3);
replaced_fix10 = nanmean(replaced_fix10,3);
replaced_time = nanmean(replaced_time,3);
replaced_time3 = nanmean(replaced_time3,3);
%%
figure
plot(moved_fixations(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of All Fixations')
ylabel('First Fixation in ROI from BCRW')
title('Moved')

figure
plot(moved_fix10(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of 1st 10 Fixations')
ylabel('First Fixation in ROI from BCRW')
title('Moved')

figure
plot(moved_time(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of All Time')
ylabel('First Fixation in ROI from BCRW')
title('Moved')

figure
plot(moved_time3(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of 1st 3 secs')
ylabel('First Fixation in ROI from BCRW')
title('Moved')

figure
plot(moved_time15(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of 0.5-1.5 secs')
ylabel('First Fixation in ROI from BCRW')
title('Moved')
%%
figure
plot(moved_fixations(1:end),moved_total(1:end),'k.')
xlabel('Proportion of All Fixations')
ylabel('Total Fixations in ROI from BCRW')
title('Moved')

figure
plot(moved_fix10(1:end),moved_total(1:end),'k.')
xlabel('Proportion of 1st 10 Fixations')
ylabel('Total Fixations in ROI from BCRW')
title('Moved')

figure
plot(moved_time(1:end),moved_total(1:end),'k.')
xlabel('Proportion of All Time')
ylabel('Total Fixations in ROI from BCRW')
title('Moved')

figure
plot(moved_time3(1:end),moved_total(1:end),'k.')
xlabel('Proportion of 1st 3 secs')
ylabel('Total Fixations in ROI from BCRW')
title('Moved')

figure
plot(moved_time15(1:end),moved_total(1:end),'k.')
xlabel('Proportion of 0.5-1.5 secs')
ylabel('Total Fixations in ROI from BCRW')
title('Moved')
%%
figure
plot(replaced_fixations(1:end),replaced_difficulty(1:end),'k.')
xlabel('Proportion of All Fixations')
ylabel('First Fixation in ROI from BCRW')
title('replaced')

figure
plot(replaced_fix10(1:end),replaced_difficulty(1:end),'k.')
xlabel('Proportion of 1st 10 Fixations')
ylabel('First Fixation in ROI from BCRW')
title('replaced')

figure
plot(replaced_time(1:end),replaced_difficulty(1:end),'k.')
xlabel('Proportion of All Time')
ylabel('First Fixation in ROI from BCRW')
title('replaced')

figure
plot(replaced_time3(1:end),replaced_difficulty(1:end),'k.')
xlabel('Proportion of 1st 3 secs')
ylabel('First Fixation in ROI from BCRW')
title('replaced')
%%
figure
plot(replaced_fixations(1:end),replaced_total(1:end),'k.')
xlabel('Proportion of All Fixations')
ylabel('Total Fixations in ROI from BCRW')
title('replaced')

figure
plot(replaced_fix10(1:end),replaced_total(1:end),'k.')
xlabel('Proportion of 1st 10 Fixations')
ylabel('Total Fixations in ROI from BCRW')
title('replaced')

figure
plot(replaced_time(1:end),replaced_total(1:end),'k.')
xlabel('Proportion of All Time')
ylabel('Total Fixations in ROI from BCRW')
title('replaced')

figure
plot(replaced_time3(1:end),replaced_total(1:end),'k.')
xlabel('Proportion of 1st 3 secs')
ylabel('Total Fixations in ROI from BCRW')
title('replaced')
