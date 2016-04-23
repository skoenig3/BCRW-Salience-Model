%Determine if various IOR_tau values result in differenece in viewing behavior
%search [#] to get to the next section
% [9] Calculate goodness of fit of BCRW for fixation location using KL-Divergence
% [10] Combine BCRW and Calculate Goodness of Fit for Fixation Location-AUC ROC
% [11] Combine BCRW and Calculate Goodness of Fit-AUC ROC Fixation Order
% [12] Calculate Salience @ BCRW Fixation Locations & Decay Rate of Salience at fixation locations
%% [9] Calculate goodness of fit of BCRW for fixation location using KL-Divergence
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
binsize = 25; %24 pixels per dva but using 25 cuz divides nicely into 600 & 800
f = fspecial('gaussian',[256,256],24);
KLnorm = NaN(36*length(image_sets),3);

IOR_dir = 'BCRW IOR TAU ';
IOR_taus = [0 50 35 25 17 12 7 3 1];

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
    
    for IOR = 1:length(IOR_taus);
        for i = 1:36;
            allfixations = zeros(imageY,imageX);
            allBCRW = zeros(imageY,imageX);
            index = (imset-1)*36+i;
            
            for t = 1:length(tags)
                if eyedatafiles(t) ~= 0;
                    
                    load([IOR_dir num2str(IOR_taus(IOR)) '\' tags{t} '-' num2str(i) '-BCRW.mat'])
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
            BCRWorder = bin2(allBCRW,binsize,binsize);
            
            binfixations(binfixations == 0) = eps;
            BCRWorder(BCRWorder == 0) = eps;
            
            binfixations = binfixations/sum(sum(binfixations));
            BCRWorder = BCRWorder/sum(sum(BCRWorder));
            
            KLnorm(index,IOR) = sum(sum(log2(binfixations./BCRWorder).*binfixations))...
                +sum(sum(log2(BCRWorder./binfixations).*BCRWorder));
        end
    end
end

pvalues = NaN(9,9);
for i = 1:9
    for ii = 1:9
        if i < ii
        [~,p] = ttest2(KLnorm(:,i),KLnorm(:,ii));
        pvalues(i,ii) = p;
        end
    end
end

allKLnorm = [];
iorindex = [];
for i = 1:9;
   allKLnorm = [allKLnorm; KLnorm(:,i)];
   iorindex = [iorindex; i*ones(size(KLnorm(:,i),1),1)];  
end
anovap = anovan(allKLnorm,iorindex);
save('IOR_KL_analysis')

figure
hold on
bar(nanmean(KLnorm))
errorb(mean(KLnorm),std(KLnorm)./sqrt(size(KLnorm,1)))
hold off
set(gca,'XTick',[1:9])
set(gca,'XTickLabel',{'IOR tau 0','IOR tau 1/50','IOR tau 1/35','IOR tau 1/25','IOR tau 1/17',...
    'IOR tau 1/12','IOR tau 1/7','IOR tau 1/3','IOR tau 1'})
ylabel('KL Divergence')
%% [10] Combine BCRW and Calculate Goodness of Fit for Fixation Location-AUC ROC
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
f = fspecial('gaussian',[256,256],24);
ROC = cell(1,length(image_sets));

IOR_dir = 'BCRW IOR TAU ';
IOR_taus = [0 50 35 25 17 12 7 3 1];

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
    for IOR = 1:length(IOR_taus)
        for i = 1:36;
            
            BCRW_at_fixations = [];
            BCRW_at_random = [];
            allBCRW = zeros(imageY,imageX);
            index = (imset-1)*36+i;
            
            for t = 1:length(tags)
                if eyedatafiles(t) ~= 0;
                    load([IOR_dir num2str(IOR_taus(IOR)) '\' tags{t} '-' num2str(i) '-BCRW.mat'])
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
                        fixBCRW = NaN(1,size(fixations,2));
                        sfixBCRW = NaN(1,size(fixations,2));
                        for iii = 1:size(fixations,2)
                            xxyy = round(fixations(:,iii));
                            xxyy(2) = imageY-xxyy(2);
                            xxyy(xxyy < 1) = 1;
                            xxyy(2,(xxyy(2) > imageY)) = imageY;
                            xxyy(1,(xxyy(1) > imageX)) = imageX;
                            fixBCRW(iii) = allBCRW(xxyy(2),xxyy(1));
                            ry = randi(600); ry(ry < 1) = 1;
                            rx = randi(800); rx(rx < 1) = 1;
                            sfixBCRW(iii) = allBCRW(ry,rx);
                        end
                        BCRW_at_fixations = [BCRW_at_fixations fixBCRW];
                        BCRW_at_random = [BCRW_at_random sfixBCRW];
                    end
                end
            end
            len = length(BCRW_at_fixations);
            thresh = 0:0.01:1;
            TP = NaN(1,length(thresh)); %True positive
            FA = NaN(1,length(thresh)); %False alarm
            for ii = 1:length(thresh)
                TP(1,ii) = sum(BCRW_at_fixations > thresh(ii))/len;
                FA(1,ii) = sum(BCRW_at_random > thresh(ii))/len;
            end
            ROC{imset}(IOR,i) = -trapz(FA(1,:),TP(1,:));
        end
    end
end

combinedROC = [];
for imset = 1:length(image_sets);
    combinedROC = [combinedROC ROC{imset}];
end

figure
hold on
bar([1:size(combinedROC,1)],mean(combinedROC,2));
errorbar(mean(combinedROC,2),std(combinedROC')/sqrt(size(combinedROC,2)),...
    '+k','linewidth',3);
set(gca,'XTick',[1:9])
set(gca,'XTickLabel',{'IOR tau 0','IOR tau 1/50','IOR tau 1/35','IOR tau 1/25','IOR tau 1/17',...
    'IOR tau 1/12','IOR tau 1/7','IOR tau 1/3','IOR tau 1'})
ylabel('AUC of ROC curve (a.u.)')
ylim([0.65 0.75])

pvalues = NaN(9,9);
for i = 1:9
    for ii = 1:9
        if i < ii
        [~,p] = ttest2(combinedROC(i,:),combinedROC(ii,:));
        pvalues(i,ii) = p;
        end
    end
end

allROC = [];
iorindex = [];
for i = 1:9;
   allROC = [allROC; combinedROC(i,:)'];
   iorindex = [iorindex; i*ones(size(combinedROC(i,:),2),1)];  
end
anovap = anovan(allROC,iorindex);
save('IOR_AUROC_analysis')
%% [11] Combine BCRW and Calculate Goodness of Fit-AUC ROC Fixation Order
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
binsize = 25; %code is optimized for binsize of 25 other values may produce fatal errors
f = fspecial('gaussian',[256,256],24);
parkhurst = {'off'}; %distamce bias for fixation order in Parkhurst, Law, Niebur (2002)
pk = 1;%index of parkhurst to use
ROC = cell(1,length(image_sets));

IOR_dir = 'BCRW IOR TAU ';
IOR_taus = [0 50 35 25 17 12 7 3 1];


BCRW_at_fixations = [];
BCRW_at_random = [];

BCRW_order_diff = cell(1,length(image_sets));
ROC = cell(length(image_sets),length(IOR_taus));
for imset = 1:length(image_sets);
    disp(['Image set-' num2str(image_sets{imset})])
    dirName = [scm_image_dir image_sets{imset}];
    cd(dirName)
    matfiles = what;
    
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
    
    for IOR = 1:length(IOR_taus)
        iBCRW_at_fixations = [];
        iBCRW_at_random = [];
        for index = 1:36;
            disp(['Running on image #' num2str(index) ' from ' image_sets{imset}])
            
            maxfixations =0;
            for t = 1:length(tags)
                load([IOR_dir num2str(IOR_taus(IOR)) '\' tags{t} '-' num2str(index) '-BCRW.mat'])
                maxfixations = max(maxfixations,max(sum(fixationtimes(:,:,1) > 0,2)));
            end
            fixorderpdf = zeros(imageY,imageX,maxfixations);
            allBCRW = zeros(imageY,imageX);
            for t = 1:length(tags)
                load([IOR_dir num2str(IOR_taus(IOR)) '\' tags{t} '-' num2str(index) '-BCRW.mat'])
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
                [i j k] = ind2sub(size(binfixorderpdf),find(binfixorderpdf == max(max(max(binfixorderpdf)))));
                if length(k) > 1
                    mind = find(k == min(k));
                    i = i(mind);
                    j = j(mind);
                    k = k(mind);
                    if length(i) > 1
                        ind = sub2ind(size(binallBCRW),i,j);
                        [~,ind2] = sort(binallBCRW(ind));
                        ind = ind(ind2(1));
                        [i j] = ind2sub([imageY/binsize,imageX/binsize],ind);
                        k = k(1);
                    end
                end
                binfixorderpdf(:,:,k) = 0;
                binfixorderpdf(i,j,:) = 0;
                BCRWorder(i,j) = k;
            end
            binallBCRW(~isnan(BCRWorder)) = NaN;
            tempBCRW = binallBCRW;
            [rr cc] = meshgrid([11 12 13 14],[15 16 17 18]);
            tempBCRW(rr,cc) = NaN;
            [rr cc] = meshgrid(1:imageX/binsize,1:imageY/binsize);
            fixnum = maxfixations+1;
            while any(any(binallBCRW > 0));
                if isnan(min(min(tempBCRW)));
                    binallBCRW(isnan(BCRWorder)) = 0;
                    BCRWorder(isnan(BCRWorder)) = fixnum;
                else
                    [i j] = find(tempBCRW == max(max(tempBCRW)));
                    if length(i) > 1;
                        ind = sub2ind(size(tempBCRW),i,j);
                        [~,ind2] = sort(binallBCRW(ind));
                        ind = ind(ind2(1));
                        [i j] = ind2sub(size(tempBCRW),ind);
                    end
                    BCRWorder(i,j) = fixnum;
                    binallBCRW(i,j) = NaN;
                    tempBCRW = binallBCRW;
                    tempBCRW(i,j) = 0;
                    if strcmpi(parkhurst{pk},'on')
                        tempBCRW = tempBCRW.*exp(-((rr-i).^2+(cc-j).^2)/(2*5^2));
                    elseif strcmpi(parkhurst{pk},'rand')
                        tempBCRW = tempBCRW.*exp(-((rr-i).^2+(cc-j).^2)/(2*randi([1 10])^2));
                    elseif strcmpi(parkhurst{pk},'off')
                        C = sqrt((rr-j).^2+(cc-i).^2)<=2;
                        tempBCRW(C == 1) = NaN;
                    end
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
                        
                        fixBCRW = NaN(1,100);
                        sfixBCRW = NaN(1,100);
                        for iii = 1:size(fixations,2)
                            xxyy = fixations(:,iii);
                            xxyy(2) = imageY-xxyy(2);
                            xxyy = round((xxyy-binsize/2)/binsize+1);
                            xxyy(xxyy < 1) = 1;
                            xxyy(2,(xxyy(2) > size(BCRWorder,1))) = size(BCRWorder,1);
                            xxyy(1,(xxyy(1) > size(BCRWorder,2))) = size(BCRWorder,2);
                            fixBCRW(iii) = abs(BCRWorder(xxyy(2),xxyy(1))-iii);
                            sfixBCRW(iii) = abs(BCRWorder(randi(numel(BCRWorder)))-iii);
                        end
                        iBCRW_at_fixations = [iBCRW_at_fixations;fixBCRW];
                        iBCRW_at_random = [iBCRW_at_random;sfixBCRW];
                    end
                end
            end
        end
        
        BCRW_at_fixations = [BCRW_at_fixations;iBCRW_at_fixations];
        BCRW_at_random = [BCRW_at_random;iBCRW_at_random];
        
        nans = find(isnan(nanmean(BCRW_at_fixations)));
        iBCRW_at_fixations(:,nans) = [];
        iBCRW_at_random(:,nans) = [];
        BCRW_at_fixations(:,nans) = [];
        BCRW_at_random(:,nans) = [];
        
        medianlen = median(sum(~isnan(BCRW_at_fixations),2));
        thresh = 0:1:numel(BCRWorder);
        TP = NaN(1,length(thresh)); %True positive
        FA = NaN(1,length(thresh)); %False alarm
        for fixnum = 1:size(iBCRW_at_fixations,2);
            for ii = 1:length(thresh)
                len = sum(~isnan(iBCRW_at_fixations(:,fixnum)));
                TP(1,ii) = sum(iBCRW_at_fixations(:,fixnum) < thresh(ii))/len;
                FA(1,ii) = sum(iBCRW_at_random(:,fixnum) < thresh(ii))/len;
            end
            ROC{imset}(IOR,fixnum) = trapz(FA(1,:),TP(1,:));
        end
        BCRW_order_diff{imset,IOR} = iBCRW_at_fixations;
    end
end

combinedROC = cell(1,9);
for i = 1:length(image_sets);
    for ii = 1:length(IOR_taus)
        combinedROC{ii} = [combinedROC{ii}; ROC{i}(ii,1:40)];
    end
end

figure
hold all
for i = 1:length(combinedROC)
    plot(mean(combinedROC{i}))
end
hold off
legend({'IOR tau 0','IOR tau 1/50','IOR tau 1/35','IOR tau 1/25','IOR tau 1/17',...
    'IOR tau 1/12','IOR tau 1/7','IOR tau 1/3','IOR tau 1'});
xlabel('Ordinal Fixation #')
ylabel('AUROC (a.u.)')


pvalues = NaN(9,9);
for i = 1:9
    for ii = 1:9
        if i < ii
            [~,p] = kstest2(combinedROC{i}(1:end),combinedROC{ii}(1:end));
            pvalues(i,ii) = p;
        end
    end
end

fixnum = [];
IOR = [];
AUROC = [];
for i = 1:9
   for set = 1:8 
    	fixnum = [fixnum; [1:40]'];
        IOR = [IOR; i*ones(40,1)];
        AUROC = [AUROC; combinedROC{i}(set,:)'];
   end
end
anova_p_AUROC = anovan(AUROC,[IOR fixnum],'model','interaction');

BCRW_order = cell(1,length(IOR_taus));
for set = 1:8
    for i = 1:9;
        BCRW_order{i} = [ BCRW_order{i}; nanmean(BCRW_order_diff{set,i}(:,1:40))]; %only want up to the median number of fixations
    end
end

figure
hold all
for i = 1:length(BCRW_order)
    plot(mean(BCRW_order{i}))
end
hold off
legend({'IOR tau 0','IOR tau 1/50','IOR tau 1/35','IOR tau 1/25','IOR tau 1/17',...
    'IOR tau 1/12','IOR tau 1/7','IOR tau 1/3','IOR tau 1'});
xlabel('Ordinal Fixation #')
ylabel('Error/Difference in Fixation Order')

fixnum = [];
IOR = [];
Difforder= [];
for i = 1:9
   for set = 1:8 
    	fixnum = [fixnum; [1:40]'];
        IOR = [IOR; i*ones(40,1)];
        Difforder = [Difforder; BCRW_order{i}(set,:)'];
   end
end
anova_p_fixorder = anovan(Difforder,[IOR fixnum],'model','interaction');

save('IOR_FixORder_analysis')
%% [12] Calculate Salience @ BCRW Fixation Locations & Decay Rate of Salience at fixation locations
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
binsize = 25; %24 pixels per dva but using 25 cuz divides nicely into 600 & 800
f = fspecial('gaussian',[256,256],50);

IOR_dir = 'BCRW IOR TAU ';
IOR_taus = [0 50 35 25 17 12 7 3 1];

salience_at_fixations = cell(1,length(IOR_taus));
for imset = 1:length(image_sets);
    disp(['Image set-' num2str(image_sets{imset})])
    dirName = [scm_image_dir image_sets{imset}];
    cd(dirName)
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
    
    
    for IOR = 1:length(IOR_taus);
        salience_at_fixations{IOR} = NaN(1152,50);
        count = 1;
        for i = 1:36;
            sal_at_fixations = NaN(100,50);
            
            for t = 1:length(tags)
                load([IOR_dir num2str(IOR_taus(IOR)) '\' tags{t} '-' num2str(i) '-BCRW.mat'])
                load(matfiles.mat{saliencemapfiles(i)},'fullmap');
                
                allBCRW = zeros(imageY,imageX);
                for ii = 1:size(fixationtimes,1);
                    tind = find(fixationtimes(ii,:,1) > 0);
                    if length(tind) > 50;
                        tind = tind(1:50);
                    end
                    for iii = 1:length(tind)
                        x = fixationtimes(ii,tind(iii),1);
                        y = fixationtimes(ii,tind(iii),2);
                        sal_at_fixations(ii,iii) = fullmap(y,x);
                    end
                end
                count = count+1;
                salience_at_fixations{IOR}(count,:) = mean(sal_at_fixations);
            end
        end
    end
end

% % MSE of scaled actual salience vs. salience @ BCRW fixation locations
% load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
%     'CombinedSalienceStatistics.mat'],'allstatistics','alldata')
% meansal = allstatistics.meanvalues(1,1,:);
% meansal = reshape(meansal,[1,40]);
% 
% MSE = [];
% scales = [];
% for i = 1:size(salience_at_fixations,2);
%     BCRWsal = nanmean(salience_at_fixations{i}(:,1:length(meansal)));
%     scale = mean(BCRWsal)/mean(meansal);
%     realsal = scale.*meansal;
%     MSE(i) = mean([BCRWsal-realsal].^2);
%     scales(i) = scale;
% end
% 
% [~,mind]=min(MSE);
% disp(['Best Fit IOR_tau is an IOR_tau of ' num2str(IOR_taus(mind))])
%%

filtsal = filtfilt(1/7*ones(1,7),1,meansal);
asymptote = min(filtsal);%use filtered version to estimate asymtote
m = meansal-asymptote; %remove offset so can fit to y = a*exp(b*x)
x = (0:39)';
y = m; 
f = fit(x,y','exp1')
plot(f,x,y)

%newer attempt
% filtsal = filtfilt(1/7*ones(1,7),1,meansal);
% asymptote = min(filtsal);%use filtered version to estimate asymtote
% zeroedsal = meansal-asymptote; 
% zeroedsal(zeroedsal <= 0) = NaN;
% logsal = log(zeroedsal);
% realind = find(~isnan(logsal));
% p = polyfit(realind,logsal(realind),1);

%old version
% realsal = meansal-min(meansal);
% zind = find(realsal == 0);
% realsal = realsal(1:zind-1);
% realsal = realsal(1:end-1);
% %     realsal = realsal/max(realsal);% don't really need to normalize cuz only looking for tau
% realsal = log(realsal);
% realsal = realsal(1:12);
% x = 1:length(realsal);
% p = polyfit(x,realsal,1);
% observed_tau = p(1);

%%
figure
colorset = hsv;
set(gca,'ColorOrder',colorset(end:-7:1,:))
hold all
for i = 1:length(IOR_taus)
    plot(nanmean(salience_at_fixations{i}))
end
%plot(scales(mind)*meansal(1:35),'k')
xlabel('Fixation Number')
ylabel('Normalized Salience')
labels= {'IOR Tau 0','IOR Tau 1/50','IOR Tau 1/35','IOR Tau 1/25',...
    'IOR Tau 1/17','IOR Tau 1/12','IOR Tau 1/7','IOR Tau 1/3','IOR Tau 1'};
legend(labels);
xlim([0 40])

load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'CombinedSalienceStatistics.mat'],'allstatistics','alldata')
meansal = allstatistics.meanvalues(1,1,:);
meansal = reshape(meansal,[1,40]);


pvalues = NaN(1,6);
[~,p] = kstest2(salience_at_fixations{1}(1:end),salience_at_fixations{2}(1:end));
pvalues(1) = p;
[~,p] = kstest2(salience_at_fixations{1}(1:end),salience_at_fixations{3}(1:end));
pvalues(2) = p;
[~,p] = kstest2(salience_at_fixations{2}(1:end),salience_at_fixations{3}(1:end));
pvalues(3) = p;

salvals = alldata(:,1,1,:);
salvals = salvals(1:end);
[~,p] = kstest2(salience_at_fixations{1}(1:35),salvals);
pvalues(4) = p;
[~,p] = kstest2(salience_at_fixations{2}(1:35),salvals);
pvalues(5) = p;
[~,p] = kstest2(salience_at_fixations{3}(1:35),salvals);
pvalues(6) = p;
% Visual Confirmation to Relate tau_IOR to decay rate of salience over fixation number
%%
slope = [];
for i = 1:length(IOR_taus)
    if IOR_taus(i) == 0 || IOR_taus(i) == 50
        val = nanmean(salience_at_fixations{i}(:,1:40));
    elseif IOR_taus(i) < 7;
        val = nanmean(salience_at_fixations{i}(:,1:7));
    else
        val = nanmean(salience_at_fixations{i}(:,1:IOR_taus(i)));
    end
    val = val-min(val);
    zind = find(val == 0);
    val = val(1:zind-1);
    val = val(1:end-1);
    %     val = val/max(val);% don't really need to normalize cuz only looking for tau
    val = log(val);
    x = 1:length(val);
    
    p = polyfit(x,val,1);
    slope(i) = p(1);
    
    if IOR_taus(i) == 0 || IOR_taus(i) == 50
        val = nanmean(salience_at_fixations{i}(:,1:40));
    elseif IOR_taus(i) < 7;
        val = nanmean(salience_at_fixations{i}(:,1:7));
    else
        val = nanmean(salience_at_fixations{i}(:,1:IOR_taus(i)));
    end
    val = val-min(val);
    x = 1:length(val);
    figure
    plot(val)
    hold on
    plot(x,exp(p(1)*x+p(2)),'r');
    if IOR_taus(i) == 0
        plot(x,exp(-3.71/50*x+p(2)),'g')
    else
        plot(x,exp(-3.71/IOR_taus(i)*x+p(2)),'g')
    end
    title(['tau_{IOR} = 1/' num2str(IOR_taus(i))])
end
%%
tau_tau = polyfit(1./IOR_taus(2:end-1),slope(2:end-1),1);
y = polyval(tau_tau,1./IOR_taus(2:end-1));
figure
hold on
plot(1./IOR_taus(2:end-1),slope(2:end-1),'k*')
plot(1./IOR_taus(2:end-1),y,'r')
hold off
xlabel('tau_{ior}')
ylabel('tau')

r2 = corrcoef(1./IOR_taus(2:end-1),slope(2:end-1));
r2 = r2(2).^2;
%%
save('Salience_at_fixation_locations')
