% Written by Seth Koenig last on 8/27/13. Re-written 9/3/15. This code
% determines how well a central bias can predict viewing behavior.

%% Calculate the average salience, image intensity, and BCRW maps
%code taken from AverageMaps.m
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
f = fspecial('gaussian',[256,256],24);

allavgsalience = zeros(imageY,imageX);
allavgfixations = zeros(imageY,imageX);
imageintensities = zeros(imageY,imageX);
allavgBCRW = zeros(imageY,imageX);
for imset = 1:length(image_sets);
    dirName = [scm_image_dir image_sets{imset}];
    disp(['Image set-' num2str(image_sets{imset})])
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

    for i = 1:36

        load([num2str(i) '-saliencemap.mat'])
        if ~any(any(isnan(fullmap)))
            allavgsalience = allavgsalience + fullmap;
        end

        img = double(rgb2gray(imread([num2str(i) '.bmp'])))+1; %values from 0-255 now 1-256
        imageintensities = imageintensities+img;
        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;

                load(['BCRW IOR TAU 17\' tags{t} '-' num2str(i) '-BCRW.mat'],'fixations')
                allavgBCRW = allavgBCRW+fixations;


                load(matfiles.mat{eyedatafiles(t)})
                fixations = fixationstats{i*2-1}.fixations;
                fixationtimes = fixationstats{i*2-1}.fixationtimes;
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
                        allavgfixations(xxyy(2),xxyy(1)) = allavgfixations(xxyy(2),xxyy(1)) + 1;
                    end
                end
            end
        end
    end
end
allavgfixations = imfilter(allavgfixations,f);
allavgBCRW = imfilter(allavgBCRW,f);
%% Using KL divergence determine the goodness of fit of the central bias
binsize = 25;

avgsalmap = bin2(allavgsalience,binsize,binsize);
avgBCRWmap = bin2(allavgBCRW,binsize,binsize);
avgfixmap = bin2(allavgfixations,binsize,binsize);
avgsalmap = avgsalmap/sum(sum(avgsalmap));
avgBCRWmap = avgBCRWmap/sum(sum(avgBCRWmap));
avgfixmap = avgfixmap/sum(sum(avgfixmap));

KLbias = NaN(288,3);

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
        index = (imset-1)*36+i;

        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;

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
        binfixations = bin2(allfixations,binsize,binsize);
        binfixations(binfixations == 0) = eps;
        binfixations = binfixations/sum(sum(binfixations));

        KLbias(index,1) = sum(sum(log2(binfixations./avgsalmap).*binfixations))...
            +sum(sum(log2(avgsalmap./binfixations).*avgsalmap));
        KLbias(index,2) = sum(sum(log2(binfixations./avgBCRWmap).*binfixations))...
            +sum(sum(log2(avgBCRWmap./binfixations).*avgBCRWmap));
        KLbias(index,3) = sum(sum(log2(binfixations./avgfixmap).*binfixations))...
            +sum(sum(log2(avgfixmap./binfixations).*avgfixmap));
    end
end

load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'Combined-KL-DiveregenceTest-CorrectedImgI.mat'],'allKL');%load previously analyzed KL divergence data
%allKL: column 1 salience vs fixation, column 2 BCRW vs fixation, column 3
%image intensity vs fixation & row by image number

%distrubations looks normal excpet a few terrible fits (very far outliers)
%so can probably assume a normal distribution and use a t-test.
KLpvalues = NaN(1,4);
[~,p] = kstest2(allKL(:,1),KLbias(:,1)); %avg salience vs salience
KLpvalues(1) = p;
[~,p] = kstest2(allKL(:,2),KLbias(:,2)); %avg BCRW vs BCRW
KLpvalues(2) = p;
[~,p] = kstest2(allKL(:,1),KLbias(:,3)); %salience vs avg fixation
KLpvalues(3) = p;
[~,p] = kstest2(allKL(:,2),KLbias(:,3)); %BCRW vs avg fixation
KLpvalues(4) = p;

pos = 2:3:8;
figure
hold on
b(1) = bar(pos(1:2)-0.5,nanmean(allKL(:,1:2)),'b');
set(b(1),'BarWidth',0.3)
b(2) = bar(pos+0.5,nanmean(KLbias),'r');
set(b(2),'BarWidth',0.3)
errorbar(pos(1:2)-0.5,nanmean(allKL(:,1:2)),nanstd(allKL(:,1:2))/sqrt(size(allKL,1)),'xk','linewidth',2);
errorbar(pos+0.5,nanmean(KLbias),nanstd(KLbias)/sqrt(size(KLbias,1)),'xk','linewidth',2);

for i = 1:length(KLpvalues);
    if i <= 2
        h = mean(KLbias(:,i)) + .25;
        pos2 = pos(i);
        if KLpvalues(i) < 0.001
            text(pos2,h,['p = ' num2str(KLpvalues(i),'%1.1e\n')],'HorizontalAlignment','center')
        else
            text(pos2,h,['p = ' num2str(KLpvalues(i),'%.3f')],'HorizontalAlignment','center')
        end
    else
        if i == 3;
            plot([2 8],[6 6],'k')
            text(5,6.05,['p = ' num2str(KLpvalues(i),'%1.1e\n')],'HorizontalAlignment','center')
        elseif i == 4;
            plot([5 8],[6.25 6.25],'k')
            text(6.5,6.3,['p = ' num2str(KLpvalues(i),'%1.1e\n')],'HorizontalAlignment','center')
        end
    end
end
hold off
legend(b,{'Normal','Central Bias'},'location','NorthEastOutside')
ylim([0 6.35])
ylabel('KL Divergence: D_{KL} (Distance in Bits)')
set(gca,'XTick',pos)
set(gca,'XTickLabel',{'Fixation vs Salience','Fixation vs BCRW',...
    'Fixation vs Average Fixation map'})

save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\Central-Bias-KL-DiveregenceTest.mat'],'KLpvalues','allKL',...
    'KLbias')
%%
%---[10] Combine BCRW and Calculate Goodness of Fit for Fixation Location-AUC ROC---%

%normalize maps to be from 0 to 1.
allavgsalience = allavgsalience-min(min(allavgsalience));
allavgsalience = allavgsalience/max(max(allavgsalience));
allavgfixations = allavgfixations-min(min(allavgfixations));
allavgfixations = allavgfixations/max(max(allavgfixations));
allavgBCRW = allavgBCRW-min(min(allavgBCRW));
allavgBCRW = allavgBCRW/max(max(allavgBCRW));

ROC = cell(1,length(image_sets));

figure
hold on
for imset = 1:length(image_sets);
    avgROC{imset} = NaN(3,36);
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
        avgsalience_at_fixations = [];
        avgsalience_at_random = [];
        avgBCRW_at_fixations = [];
        avgBCRW_at_random = [];
        avgfix_at_fixations = [];
        avgfix_at_random = [];
        
        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;
                
                load(matfiles.mat{eyedatafiles(t)})
                fixations = fixationstats{i*2-1}.fixations;
                if ~isempty(fixations)
                    if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                            fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                        fixations(:,1) = [];
                    end
                    
                    
                    avgfixsal  = NaN(1,size(fixations,2));
                    savgfixsal = NaN(1,size(fixations,2));
                    avgfixBCRW = NaN(1,size(fixations,2));
                    savgfixBCRW = NaN(1,size(fixations,2));
                    avgfixfix  = NaN(1,size(fixations,2));
                    savgfixfix = NaN(1,size(fixations,2));
                    
                    for iii = 1:size(fixations,2)
                        xxyy = round(fixations(:,iii));
                        xxyy(2) = imageY-xxyy(2);
                        xxyy(xxyy < 1) = 1;
                        xxyy(2,(xxyy(2) > imageY)) = imageY;
                        xxyy(1,(xxyy(1) > imageX)) = imageX;
                        
                        avgfixsal(iii) = allavgsalience(xxyy(2),xxyy(1));
                        avgfixBCRW(iii) = allavgBCRW(xxyy(2),xxyy(1));
                        avgfixfix(iii) = allavgfixations(xxyy(2),xxyy(1));
                        
                        ry = randi(600); ry(ry < 1) = 1;
                        rx = randi(800); rx(rx < 1) = 1;
                        
                        savgfixsal(iii) = allavgsalience(ry,rx);
                        savgfixBCRW(iii) = allavgBCRW(ry,rx);
                        savgfixfix(iii) = allavgfixations(ry,rx);
                    end
                    avgsalience_at_fixations = [avgsalience_at_fixations avgfixsal];
                    avgsalience_at_random = [avgsalience_at_random savgfixsal];
                    avgBCRW_at_fixations = [avgBCRW_at_fixations avgfixBCRW];
                    avgBCRW_at_random = [avgBCRW_at_random savgfixBCRW];
                    avgfix_at_fixations = [avgfix_at_fixations avgfixfix];
                    avgfix_at_random = [avgfix_at_random savgfixfix];
                end
            end
        end
        
        len = length(avgsalience_at_fixations);
        thresh = 0:0.01:1;
        TP = NaN(3,length(thresh)); %True positive
        FA = NaN(3,length(thresh)); %False alarm
        for ii = 1:length(thresh)
            TP(1,ii) = sum(avgsalience_at_fixations > thresh(ii))/len;
            FA(1,ii) = sum(avgsalience_at_random > thresh(ii))/len;
            TP(2,ii) = sum(avgBCRW_at_fixations > thresh(ii))/len;
            FA(2,ii) = sum(avgBCRW_at_random > thresh(ii))/len;
            TP(3,ii) = sum(avgfix_at_fixations > thresh(ii))/len;
            FA(3,ii) = sum(avgfix_at_random > thresh(ii))/len;
        end
        avgROC{imset}(1,i) = -trapz(FA(1,:),TP(1,:));
        avgROC{imset}(2,i) = -trapz(FA(2,:),TP(2,:));
        avgROC{imset}(3,i) = -trapz(FA(3,:),TP(3,:));
        plot(FA(1,:),TP(1,:),'r')
        plot(FA(2,:),TP(2,:),'g')
        plot(FA(3,:),TP(3,:),'b')
    end
end
plot(thresh,thresh,'k--','linewidth',3)
hold off
legend('Avg Salience','Avg BCRW','Avg Fixation','location','Northeastoutside')
xlabel('False Alarm Rate (FP)')
ylabel('True Positive (TP)')

load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'Combined-ROC-Corrected-ImageI.mat'],'ROC')

combinedROC = [];
avgcombinedROC = [];
for imset = 1:length(image_sets);
    combinedROC = [combinedROC ROC{imset}]; 
    avgcombinedROC = [avgcombinedROC avgROC{imset}];
end


pos = 2:3:8;
figure
hold on
b(1) = bar(pos(1:2)-0.5,nanmean(combinedROC(1:2,:),2),'b');
set(b(1),'BarWidth',0.3)
b(2) = bar(pos+0.5,nanmean(avgcombinedROC'),'r');
set(b(2),'BarWidth',0.3)
errorbar(pos(1:2)-0.5,nanmean(combinedROC(1:2,:),2),nanstd(combinedROC(1:2,:)')/sqrt(size(combinedROC,2)),'xk','linewidth',2);
errorbar(pos+0.5,nanmean(avgcombinedROC'),nanstd(avgcombinedROC')/sqrt(size(avgcombinedROC',1)),'xk','linewidth',2);
legend(b,{'Normal','Central Bias'},'location','NorthEastOutside')
ylim([0.5 0.75])
ylabel('AUROC (a.u.)')
set(gca,'XTick',pos)
set(gca,'XTickLabel',{'Fixation vs Salience','Fixation vs BCRW',...
    'Fixation vs Average Fixation map'})

zpvalues = NaN(1,size(avgcombinedROC,1));
for z = 1:3
    [~,p] = ztest(avgcombinedROC(z,:),0.5,std(avgcombinedROC(z,:)));
    zpvalues(z) = p;
end

kspvalues = NaN(1,size(combinedROC,1)+1); %not all distributions look normal especially avg fixation map
[~,p] = kstest2(combinedROC(1,:),avgcombinedROC(1,:));
kspvalues(1) = p;
[~,p] = kstest2(combinedROC(2,:),avgcombinedROC(2,:));
kspvalues(2) = p;
[~,p] = kstest2(combinedROC(1,:),avgcombinedROC(3,:));
kspvalues(3) = p;
[~,p] = kstest2(combinedROC(2,:),avgcombinedROC(3,:));
kspvalues(4) = p;

for i = 1:3
    text(pos(i)+0.5,mean(avgcombinedROC(i,:))+0.02,['p = ' num2str(p,'%1.1e\n')],'HorizontalAlignment','center')
end
for i = 1:length(kspvalues);
    if i <= 2
        h = mean(combinedROC(i,:)) + .02;
        pos2 = pos(i);
        if kspvalues(i) < 0.001
            text(pos2,h,['p = ' num2str(kspvalues(i),'%1.1e\n')],'HorizontalAlignment','center')
        else
            text(pos2,h,['p = ' num2str(kspvalues(i),'%.3f')],'HorizontalAlignment','center')
        end
    else
        if i == 3;
            plot([2 8],[0.75 0.75],'k')
            text(5,0.76,['p = ' num2str(kspvalues(i),'%1.1e\n')],'HorizontalAlignment','center')
        elseif i == 4;
            plot([5 8],[0.79 0.79],'k')
            text(6.5,0.8,['p = ' num2str(kspvalues(i),'%1.1e\n')],'HorizontalAlignment','center')
        end
    end
end
hold off
ylim([0.5 0.825])
%%
save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\'...
    'SCM Image Sets\Central-Bias-AUROC.mat'],'avgROC','zpvalues','kspvalues')