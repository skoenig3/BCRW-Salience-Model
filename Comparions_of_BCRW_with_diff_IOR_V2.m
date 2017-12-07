%updated version written 5/15/2017 SDK, meant to be streamlined and more
%efficeint so easier to understand and faster to run
%% [9] Calculate goodness of fit of BCRW for fixation location using KL-Divergence
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
binsize = 25; %24 pixels per dva but using 25 cuz divides nicely into 600 & 800
f = fspecial('gaussian',[155,155],24); %1 dva smoothing

IOR_dir = 'BCRW IOR TAU Simulations\';
IOR_taus = [0 35 30 25 20 17 14 12 9 7 5 3 1];

KLnorm = NaN(36*length(image_sets),length(IOR_taus));
ROC = NaN(36*length(image_sets),length(IOR_taus));

salience_at_fixations = cell(1,length(IOR_taus));
median_BCRW_fixation_count =  cell(1,length(IOR_taus));
for it = 1:length(IOR_taus)
    salience_at_fixations{it} = NaN(36*6,50);
    median_BCRW_fixation_count{it} = NaN(1,36*6);
end

for imset = 1:length(image_sets);
    disp(['Image set-' num2str(image_sets{imset})])
    dirName = [scm_image_dir image_sets{imset} '\'];
    
    matfiles = what(dirName);
    
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
            load([dirName num2str(i) '-saliencemap.mat'],'fullmap');
            allfixations = zeros(imageY,imageX);
            allBCRW = zeros(imageY,imageX);
            BCRW_at_fixations = [];
            BCRW_at_random = [];
            index = (imset-1)*36+i;
            
            sal_at_fixations = NaN(400,50);
            count = 1;
            for t = 1:length(tags)
                if IOR_taus(IOR) == 0
                    load([dirName IOR_dir tags{t} '-' num2str(i) '-' num2str(IOR_taus(IOR)) '-BCRW.mat'])
                else
                    load([dirName IOR_dir tags{t} '-' num2str(i) '-' num2str(1/IOR_taus(IOR)) '-BCRW.mat'])
                end
                allBCRW = allBCRW+fixations;
                BCRWfixations = fixationtimes;
                
                %---Calculate Salience at BCRW Fixations---%
                for ii = 1:size(BCRWfixations,1);
                    tind = find(BCRWfixations(ii,:,1) > 0);
                    if length(tind) > 50;
                        tind = tind(1:50);
                    end
                    for iii = 1:length(tind)
                        x = BCRWfixations(ii,tind(iii),1);
                        y = BCRWfixations(ii,tind(iii),2);
                        sal_at_fixations(count,iii) = fullmap(y,x);
                    end
                    count = count+1;
                end
                
                %---Calculate PDF of observed and Predicted Fixations---%
                load([dirName matfiles.mat{eyedatafiles(t)}])
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
            
            allBCRW_ROC = allBCRW;
            
            %---for KL Divergence analysis---%
            %smooth
            allfixations = imfilter(allfixations,f);
            allBCRW = imfilter(allBCRW,f);
            
            %bin
            binfixations = bin2(allfixations,binsize,binsize);
            binBCRW = bin2(allBCRW,binsize,binsize);
            
            %replace zeros with eps
            binfixations(binfixations == 0) = eps;
            binBCRW(binBCRW == 0) = eps;
            
            %turn to PDFS
            binfixations = binfixations/sum(sum(binfixations));
            binBCRW = binBCRW/sum(sum(binBCRW));
            
            KLnorm(index,IOR) = sum(sum(log2(binfixations./binBCRW).*binfixations))...
                +sum(sum(log2(binBCRW./binfixations).*binBCRW));
            
            %---for ROC analysis---%
            
            %normalize 0 to 1  for ROC analysis
            allBCRW_ROC = imfilter(allBCRW_ROC,f);
            allBCRW_ROC = allBCRW_ROC - min(allBCRW_ROC(:));
            allBCRW_ROC = allBCRW_ROC/max(allBCRW_ROC(:));
            
            for t = 1:length(tags)
                if eyedatafiles(t) ~= 0;
                    load([dirName matfiles.mat{eyedatafiles(t)}])
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
            ROC(index,IOR) = -trapz(FA(1,:),TP(1,:));
            
            salience_at_fixations{IOR}(index,:) = nanmean(sal_at_fixations);
            median_BCRW_fixation_count{IOR}(index) = median(sum(~isnan(sal_at_fixations')));
        end
    end
end
%%
%---Statistical Analysis of KL Divergence---%
pvalues = NaN(length(IOR_taus),length(IOR_taus));
xlabl = {};
for i = 1:length(IOR_taus)
    for ii = 1:length(IOR_taus)
        if i < ii
            [~,p] = ttest2(KLnorm(:,i),KLnorm(:,ii));
            pvalues(i,ii) = p;
        end
    end
   xlabl = [xlabl {['tau 1/' num2str(IOR_taus(i))]}];
end

allKLnorm = [];
iorindex = [];
for i = 1:length(IOR_taus)
    allKLnorm = [allKLnorm; KLnorm(:,i)];
    iorindex = [iorindex; i*ones(size(KLnorm(:,i),1),1)];
end
[P,T,STATS,TERMS] = anovan(allKLnorm,iorindex);
COMPARISON = multcompare(STATS);
pkw = kruskalwallis(allKLnorm,iorindex);

figure
hold on
bar(nanmean(KLnorm))
errorb(mean(KLnorm),std(KLnorm)./sqrt(size(KLnorm,1)))
hold off
set(gca,'XTick',[1:length(IOR_taus)])
set(gca,'XTickLabel',xlabl)
ylabel('KL Divergence')

%%
%---Statistical Analysis of ROC---%
figure
hold on
bar([1:size(ROC,2)],mean(ROC,1));
errorbar(mean(ROC,1),std(ROC)/sqrt(size(ROC,1)),...
    '+k','linewidth',3);
ylabel('AUC of ROC curve (a.u.)')
ylim([0.65 0.75])
set(gca,'XTick',[1:length(IOR_taus)])
set(gca,'XTickLabel',xlabl)

pvalues = NaN(length(IOR_taus),length(IOR_taus));
for i = 1:length(IOR_taus)
    for ii = 1:length(IOR_taus)
        if i < ii
            [~,p] = ttest2(ROC(i,:),ROC(ii,:));
            pvalues(i,ii) = p;
        end
    end
end

figure
allROC = [];
iorindex = [];
for i = 1:length(IOR_taus);
    allROC = [allROC; ROC(:,i)];
    iorindex = [iorindex; i*ones(size(ROC,1),1)];
end
[P,T,STATS,TERMS] = anovan(allROC,iorindex);
COMPARISON = multcompare(STATS);
%%
median_num_fix = median(cell2mat(median_BCRW_fixation_count));
%---Results of Saliecne at Fixation Analysis---%
figure
hold all
for i = 1:length(IOR_taus)
    plot(nanmean(salience_at_fixations{i}))
end
xlabel('Fixation Number')
ylabel('Normalized Salience (a.u.)')
xlim([0 median_num_fix])

% load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
%     'CombinedSalienceStatistics.mat'],'allstatistics','alldata')
% meansal = allstatistics.meanvalues(1,1,:);
% meansal = reshape(meansal,[1,40]);
%%
slope = NaN(1,length(IOR_taus));
for i = 1:length(IOR_taus)
    if IOR_taus(i) <= 2 || IOR_taus(i) > median_num_fix
        continue        
    end
    val = nanmean(salience_at_fixations{i}(:,1:median_num_fix));   

    val = val-min(val);
    zind = find(val == 0);
    val = val(1:zind-1);
    val = val(1:end-1);
    %     val = val/max(val);% don't really need to normalize cuz only looking for tau
    val = log(val);
    x = 1:length(val);
    
    p = polyfit(x,val,1);
    slope(i) = p(1);
    
  
    val = val-min(val);
    x = 1:length(val);
    figure
    hold on
    plot(x,exp(p(1)*x+p(2)),'r');
    if IOR_taus(i) == 0
        plot(x,exp(-2.7/50*x+p(2)),'g')
    else
        plot(x,exp(-2.7/IOR_taus(i)*x+p(2)),'g')
    end
    title(['tau_{IOR} = 1/' num2str(IOR_taus(i))])
end

y2fit = IOR_taus(~isnan(slope));
x2fit = slope(~isnan(slope));
tau_tau = polyfit(1./y2fit,x2fit,1);
y = polyval(tau_tau,1./y2fit);
figure
hold on
plot(1./IOR_taus,slope,'k*')
plot(1./y2fit,y,'r')
hold off
xlabel('tau_{ior}')
ylabel('tau')

r2 = corrcoef(1./y2fit,x2fit);
r2 = r2(2).^2;
%%
save('BCRW_Different_IOR_analysis_Ignore_3')
%% Calculate Percentage of Image Covered by BCRW
% scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
% image_sets = {'Set006','Set007','Set008','Set009',...
%     'SetE001','SetE002','SetE003','SetE004'};
% tags = {'MP','TT','JN','IW'};
% imageX = 800;
% imageY = 600;
% nml = imageX*imageY;
% binsize = 25; %24 pixels per dva but using 25 cuz divides nicely into 600 & 800
% IOR_dir = 'BCRW IOR TAU Simulations\';
% IOR_taus = [0 35 30 25 20 17 14 12 9 7 5 3 1];
% 
% BCRW_percent_coverage = NaN(36*length(image_sets),length(IOR_taus));
% 
% IOR_area = 48;%2 dva
% [rr,cc] = meshgrid(1:imageX,1:imageY);
% 
% for imset = 1:length(image_sets);
%     for IOR = 1:length(IOR_taus);
%         dirName = [scm_image_dir image_sets{imset} '\'];
%         for i = 1:36;
%             index = (imset-1)*36+i;
%             coverage = NaN(100,length(tags));
%             for t = 1:length(tags)
%                 if IOR_taus(IOR) == 0
%                     load([dirName IOR_dir tags{t} '-' num2str(i) '-' num2str(IOR_taus(IOR)) '-BCRW.mat'])
%                 else
%                     load([dirName IOR_dir tags{t} '-' num2str(i) '-' num2str(1/IOR_taus(IOR)) '-BCRW.mat'])
%                 end
%                 disp(['Image set-' image_sets{imset} ' IOR_tau = ' num2str(IOR_taus(IOR)) ...
%                     ' Image# ' num2str(i) ' Monkey ' tags{t}])
%                 parfor ii = 1:size(fixationtimes,1);
%                     img = zeros(imageY,imageX);
%                     tind = find(fixationtimes(ii,:,1) > 0);
%                     for iii = 1:length(tind)
%                         x = fixationtimes(ii,tind(iii),1);
%                         y = fixationtimes(ii,tind(iii),2);
%                         C = sqrt((rr-x).^2+(cc-y).^2)<=IOR_area;
%                         img(C) = 1;
%                     end
%                     coverage(ii,t) = sum(sum(img))/nml;
%                 end
%             end
%             BCRW_percent_coverage(index,IOR) = mean(coverage(:));
%         end
%     end
% end
% save('BCRWCoverage')
% 
% %%
% 
% figure
% allcoverage = [];
% iorindex = [];
% for i = 1:length(IOR_taus);
%     allcoverage = [allcoverage; BCRW_percent_coverage(:,i)];
%     iorindex = [iorindex; i*ones(size(BCRW_percent_coverage,1),1)];
% end
% [P,T,STATS,TERMS] = anovan(allcoverage,iorindex);
% COMPARISON = multcompare(STATS);
% pk = kruskalwallis(allcoverage,iorindex);