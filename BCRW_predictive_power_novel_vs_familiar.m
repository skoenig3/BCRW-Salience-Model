%---[10] Combine BCRW and Calculate Goodness of Fit for Fixation Location-AUC ROC---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','JN','IW'};
imageX = 800;
imageY = 600;
f = fspecial('gaussian',[256,256],24);
ROC = cell(1,length(image_sets));
for imset = 1:length(image_sets);
    ROC{imset} = NaN(3,36,2);
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
    
    matfiles2 = what('BCRW IOR TAU 17');
    BCRWfiles = NaN(length(tags),36);
    for i = 1:length(matfiles2.mat)
        if ~isempty(strfind(matfiles2.mat{i},'BCRW'))
            for t = 1:length(tags)
                if ~isempty(strfind(matfiles2.mat{i},tags{t}))
                    dashes = strfind(matfiles2.mat{i},'-');
                    if ~isempty(str2num(matfiles2.mat{i}(dashes(2)-1)))
                        BCRWfiles(t,str2double(matfiles2.mat{i}(dashes(1)+1:dashes(2)-1))) = i;
                    end
                    break
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
        rsalience_at_fixations = [];
        rsalience_at_random = [];
        rBCRW_at_fixations = [];
        rBCRW_at_random = [];
        rI_at_fixations = [];
        rI_at_random = [];
        allBCRW = zeros(imageY,imageX);
        index = (imset-1)*36+i;
        
        load(matfiles.mat{saliencemapfiles(i)})
        allsalience = fullmap;
        
        imageintensities = 255-double(rgb2gray(imread([num2str(i) '.bmp'])));
        imageintensities = imageintensities/max(max(imageintensities));
        
        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;
                load(['BCRW IOR TAU 17\' matfiles2.mat{BCRWfiles(t,i)}],'fixations')
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
                end
                
                if length(fixationstats) >= i*2 && ~isempty(fixationstats{i*2})
                    repfixations = fixationstats{i*2}.fixations;
                    if ~isempty(repfixations)
                        if repfixations(1,1) > imageX/2-100 && repfixations(1,1) < imageX/2+100 &&...
                                repfixations(2,1) < imageY/2+100 && repfixations(2,1) > imageY/2-100
                            repfixations(:,1) = [];
                        end
                    end
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
                
                rfixsal  = NaN(1,size(repfixations,2));
                rsfixsal = NaN(1,size(repfixations,2));
                rfixBCRW = NaN(1,size(repfixations,2));
                rsfixBCRW = NaN(1,size(repfixations,2));
                rfixI = NaN(1,size(repfixations,2));
                rsfixI = NaN(1,size(repfixations,2));
                for iii = 1:size(repfixations,2)
                    xxyy = round(repfixations(:,iii));
                    xxyy(2) = imageY-xxyy(2);
                    xxyy(xxyy < 1) = 1;
                    xxyy(2,(xxyy(2) > imageY)) = imageY;
                    xxyy(1,(xxyy(1) > imageX)) = imageX;
                    rfixsal(iii) = allsalience(xxyy(2),xxyy(1));
                    rfixI(iii) = imageintensities(xxyy(2),xxyy(1));
                    rfixBCRW(iii) = allBCRW(xxyy(2),xxyy(1));
                    ry = randi(600); ry(ry < 1) = 1;
                    rx = randi(800); rx(rx < 1) = 1;
                    rsfixsal(iii) = allsalience(ry,rx);
                    rsfixI(iii) = imageintensities(ry,rx);
                    rsfixBCRW(iii) = allBCRW(ry,rx);
                end
                salience_at_fixations = [salience_at_fixations fixsal];
                salience_at_random = [salience_at_random sfixsal];
                BCRW_at_fixations = [BCRW_at_fixations fixBCRW];
                BCRW_at_random = [BCRW_at_random sfixBCRW];
                I_at_fixations = [I_at_fixations fixI];
                I_at_random = [I_at_random sfixI];
                rsalience_at_fixations = [rsalience_at_fixations rfixsal];
                rsalience_at_random = [rsalience_at_random rsfixsal];
                rBCRW_at_fixations = [rBCRW_at_fixations rfixBCRW];
                rBCRW_at_random = [rBCRW_at_random rsfixBCRW];
                rI_at_fixations = [rI_at_fixations rfixI];
                rI_at_random = [rI_at_random rsfixI];
            end
        end
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
        ROC{imset}(1,i,1) = -trapz(FA(1,:),TP(1,:));
        ROC{imset}(2,i,1) = -trapz(FA(2,:),TP(2,:));
        ROC{imset}(3,i,1) = -trapz(FA(3,:),TP(3,:));
        
        len = length(rsalience_at_fixations);
        thresh = 0:0.01:1;
        TP = NaN(3,length(thresh)); %True positive
        FA = NaN(3,length(thresh)); %False alarm
        for ii = 1:length(thresh)
            TP(1,ii) = sum(rsalience_at_fixations > thresh(ii))/len;
            FA(1,ii) = sum(rsalience_at_random > thresh(ii))/len;
            TP(2,ii) = sum(rBCRW_at_fixations > thresh(ii))/len;
            FA(2,ii) = sum(rBCRW_at_random > thresh(ii))/len;
            TP(3,ii) = sum(rI_at_fixations > thresh(ii))/len;
            FA(3,ii) = sum(rI_at_random > thresh(ii))/len;
        end
        ROC{imset}(1,i,2) = -trapz(FA(1,:),TP(1,:));
        ROC{imset}(2,i,2) = -trapz(FA(2,:),TP(2,:));
        ROC{imset}(3,i,2) = -trapz(FA(3,:),TP(3,:));
    end
end
%%
combinedROC = [];
for imset = 1:length(image_sets);
    combinedROC = [combinedROC ROC{imset}];
end
%%
novsal = combinedROC(1,:,1);
novsal = novsal(1:end);
repsal = combinedROC(1,:,2);
repsal = repsal(1:end);
[~,novsalzp] = ztest(novsal,0.5,std(novsal))
[~,repsalzp] = ztest(repsal,0.5,std(repsal))
[~,salp] = ttest2(repsal,novsal)
%%
novBCRW = combinedROC(2,:,1);
novBCRW = novBCRW(1:end);
repBCRW = combinedROC(2,:,2);
repBCRW = repBCRW(1:end);
[~,novBCRWzp] = ztest(novBCRW,0.5,std(novBCRW))
[~,repBCRWzp] = ztest(repBCRW,0.5,std(repBCRW))
[~,BCRWp] = ttest2(repBCRW,novBCRW)
%%
novI = combinedROC(3,:,1);
novI = novI(1:end);
repI = combinedROC(3,:,2);
repI = repI(1:end);
[~,novIzp] = ztest(novI,0.5,std(novI))
[~,repIzp] = ztest(repI,0.5,std(repI))
[~,Ip] = ttest2(repI,novI)
%%
figure
bar([mean(novsal) mean(repsal) mean(novBCRW) mean(repBCRW)...
    mean(novI) mean(repI)]);
hold on
errorbar([mean(novsal) mean(repsal) mean(novBCRW) mean(repBCRW)...
    mean(novI) mean(repI)],[std(novsal) std(repsal) std(novBCRW) std(repBCRW)...
    std(novI) std(repI)]/sqrt(length(repI)),'k.')
hold off
ylabel('AUROC (a.u)')