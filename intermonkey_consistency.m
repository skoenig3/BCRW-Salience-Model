%---[10] Combine BCRW and Calculate Goodness of Fit for Fixation Location-AUC ROC---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
f = fspecial('gaussian',[256,256],24);
ROC = cell(1,length(image_sets));

for imset = 1:length(image_sets);
    ROC{imset} = NaN(4,36);
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
        flag = 0;
        allfixations = cell(1,length(tags));
        fixationmaps = zeros(imageY,imageX,length(tags));
        for t = 1:length(tags)
            load(matfiles.mat{eyedatafiles(t)})
            fixations = fixationstats{i*2-1}.fixations;
            if ~isempty(fixations)
                fixationtimes = fixationstats{i*2-1}.fixationtimes;
                if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                        fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                    fixations(:,1) = [];
                end
                allfixations{t} = fixations;
                for iii = 1:length(fixations)
                    x = round(fixations(1,iii));
                    y = round(fixations(2,iii));
                    x(x < 1) = 1; x(x > imageX) = imageX;
                    y(y < 1) = 1; y(y > imageY) = imageY;
                    fixationmaps(y,x,t) =  fixationmaps(y,x,t)+1;
                end
            else
                flag = 1;
                break
            end
        end
        if flag == 1; %if one image is missing fixations can't really do this analysis
            continue
        end
        
        for tt = 1:length(tags);
            fixationmaps(:,:,tt) = imfilter(fixationmaps(:,:,tt),f);
        end
        
        %get the map for compairing all others
        allfixationmaps = zeros(imageY,imageX,length(tags));
        for tt = 1:length(tags);
            for ttt = 1:length(tags)
                if tt ~= ttt
                    allfixationmaps(:,:,tt) = allfixationmaps(:,:,tt)+fixationmaps(:,:,ttt);
                end
            end
            allfixationmaps(:,:,tt) = allfixationmaps(:,:,tt) - min(min(allfixationmaps(:,:,tt))); %normalize 0 to 1
            allfixationmaps(:,:,tt) = allfixationmaps(:,:,tt)/max(max(allfixationmaps(:,:,tt)));
        end
        
        for t = 1:length(tags)
            fixations = allfixations{t};
            for iii = 1:size(fixations,2)
                xxyy = round(fixations(:,iii));
                xxyy(xxyy < 1) = 1;
                xxyy(2,(xxyy(2) > imageY)) = imageY;
                xxyy(1,(xxyy(1) > imageX)) = imageX;
                ry = randi(600); ry(ry < 1) = 1;
                rx = randi(800); rx(rx < 1) = 1;
               
                fixp(iii) = allfixationmaps(xxyy(2),xxyy(1),t);
                sfixp(iii) = allfixationmaps(ry,rx,t);
            end
            
            len = length(fixp);
            thresh = 0:0.01:1;
            TP = NaN(4,length(thresh)); %True positive
            FA = NaN(4,length(thresh)); %False alarm
            for ii = 1:length(thresh)
                TP(1,ii) = sum(fixp > thresh(ii))/len;
                FA(1,ii) = sum(sfixp > thresh(ii))/len;
            end
            ROC{imset}(t,i) = -trapz(FA(1,:),TP(1,:));
        end
    end
end

combinedROC = [];
for i = 1:length(ROC)
    combinedROC = [combinedROC ROC{i}];
end

%%
zpvalues = NaN(1,length(tags));
for i = 1:length(tags)
    [~,p] = ztest(combinedROC(i,:),0.5,nanstd(combinedROC(i,:)),'tail','right');
    zpvalues(i) = p;
end
[~,zpall] = ztest(combinedROC(1:end),0.5,nanstd(combinedROC(1:end)));

monkey = ones(size(combinedROC));
monkey(2,:) = 2*monkey(2,:);
monkey(3,:) = 3*monkey(3,:);
monkey(4,:) = 4*monkey(4,:);
[P,ANOVATAB,STATS] = anova1(combinedROC(1:end),monkey(1:end));
multcompare(STATS);
%% Compared AUROC for intermonkey consistency to BCRW salience and Img I
load('C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\Combined-ROC-Corrected-ImageI.mat')

combinedROC2 = []; %for BCRW
for imset = 1:length(image_sets);
    combinedROC2 = [combinedROC2 ROC{imset}];
end

tpvalues = NaN(1,3);
for i = 1:3
    [~,p] = ttest2(combinedROC2(i,:),mean(combinedROC));
    tpvalues(i) = p;
end