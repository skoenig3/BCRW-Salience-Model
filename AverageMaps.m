%% Code calculates the average salience, image intensity, and BCRW maps
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
f = fspecial('gaussian',[256,256],24);

allsalience = zeros(imageY,imageX);
allfixations = zeros(imageY,imageX);
imageintensities = zeros(imageY,imageX);
allBCRW = zeros(imageY,imageX);
allCRW = zeros(imageY,imageX);
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
            allsalience = allsalience + fullmap;
        end
        
        img = double(rgb2gray(imread([num2str(i) '.bmp'])))+1; %values from 0-255 now 1-256
        img = img-min(min(img)); %zero
        img = img/max(max(img)); %scale to 1
        imageintensities = imageintensities+img;
        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;
                
                load(['BCRW IOR TAU 17\' tags{t} '-' num2str(i) '-BCRW.mat'],'fixations')
                allBCRW = allBCRW+fixations;
                                
                load(['CRW IOR TAU 17\' tags{t} '-' num2str(i) '-CRW.mat'],'fixations')
                allCRW = allCRW+fixations;
                
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
                        allfixations(xxyy(2),xxyy(1)) = allfixations(xxyy(2),xxyy(1)) + 1;
                    end
                end
            end
        end
    end
end
%%
allfixations = imfilter(allfixations,f);
allBCRW = imfilter(allBCRW,f);
allCRW = imfilter(allCRW,f);