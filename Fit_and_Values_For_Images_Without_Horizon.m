% clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';

imageX = 800;
imageY = 600;
%  
%load ROC analysis data
load([data_dir 'Combined-ROC-Corrected-ImageI.mat']);

%load KL divergence analysis data
load([data_dir 'Combined-KL-DiveregenceTest-CorrectedImgI.mat']);

threshold = 1/3*imageY;

%col 1 no horizon, col 2 with horizon  > 1/3 of image
ROC_values = cell(1,2);
KL_values = cell(1,2);
sal_values = cell(1,2);
imgI_values = cell(1,2); 
% 
% %%---[4] Calculate Average Salience at Each Fixation Across mutliple data sets---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};

minlen = 35; %median number of fixations according to viewing behavior section. 

%manually determine which images contain a lot of sky
all_horizons = cell(1,length(image_sets));
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])

    horizon = zeros(1,36);
    for img = 1:36
        imag = imread([num2str(img) '.bmp']);
        
        imshow(imag);
        hold on
        plot([0 800],[threshold threshold],'r--')
        hold off
        reply = input('Does image contain horizon below threshold [y,n]');
        if strcmpi(reply,'y') || reply == 1
            horizon(img) = 1;
        end
        
    end
    
    all_horizons{SET} = horizon;
    
end
save([data_dir 'all_horizons.mat','all_horizons','threshold'])
%% Calculate the image intesnity, salience and fixation maps for images with and without horizons

load([data_dir 'all_horizons.mat'],'all_horizons')
horizon_imageintensities = zeros(imageY,imageX);
no_horizon_imageintensities = zeros(imageY,imageX);

horizon_salmaps = zeros(imageY,imageX);
no_horizon_salmaps = zeros(imageY,imageX);

horizon_allfixations = zeros(imageY,imageX);
no_horizon_allfixations = zeros(imageY,imageX);

horizon_BCRW = zeros(imageY,imageX);
no_horizon_BCRW = zeros(imageY,imageX);
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
        
          img = double(rgb2gray(imread([num2str(i) '.bmp'])))+1; %values from 0-255 now 1-256
          img = img-min(min(img)); %zero
          img = img/max(max(img)); %scale to 1
          
          load([num2str(i) '-saliencemap.mat'])
          if any(any(isnan(fullmap)))
              fullmap = 0;
          end
          
          if all_horizons{imset}(i) == 1 %has a horizon
              horizon_imageintensities = horizon_imageintensities+img;
              horizon_salmaps = horizon_salmaps+fullmap;
          else
              no_horizon_imageintensities = no_horizon_imageintensities+img;
              no_horizon_salmaps = no_horizon_salmaps+fullmap;
          end
          
          for t = 1:length(tags)
              if eyedatafiles(t) ~= 0;
                  
                  try
                      load(['BCRW IOR TAU 17\' tags{t} '-' num2str(i) '-BCRW.mat'],'fixations')
                  catch
                      continue
                  end
                  
                  if all_horizons{imset}(i) == 1 %has a horizon
                      horizon_BCRW = horizon_BCRW+fixations;
                  else
                      no_horizon_BCRW = no_horizon_BCRW+fixations;
                  end
                  
                  
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
                          
                          if all_horizons{imset}(i) == 1 %has a horizon
                              horizon_allfixations(xxyy(2),xxyy(1)) = horizon_allfixations(xxyy(2),xxyy(1))+1;
                          else
                               no_horizon_allfixations(xxyy(2),xxyy(1)) = no_horizon_allfixations(xxyy(2),xxyy(1))+1;
                          end
                      end
                  end
              end
          end
    end
end

%%
num_horizon = sum(sum(cell2mat(all_horizons)));
num_no_horizon = sum(sum(cell2mat(all_horizons)==0));

clims = NaN(2,2);

figure
subplot(2,2,1)
imagesc(horizon_imageintensities./num_horizon)
colormap('jet')
axis off
axis equal
title('Image Intesnity: Images with Horizon')
clims(1,:) = caxis;


subplot(2,2,2)
imagesc(no_horizon_imageintensities./num_no_horizon)
colormap('jet')
axis off
axis equal
title('Image Intesnity: Images without Horizon')
clims(2,:) = caxis;

minc = min(clims(:,1));
maxc = max(clims(:,2));

for sb = 1:2
    subplot(2,2,sb)
    caxis([minc maxc])
end

clims = NaN(2,2);

subplot(2,2,3)
imagesc(horizon_salmaps./num_horizon)
colormap('jet')
axis off
axis equal
title('Salience: Images with Horizon')
clims(1,:) = caxis;

subplot(2,2,4)
imagesc(no_horizon_salmaps./num_no_horizon)
colormap('jet')
axis off
axis equal
title('Salience: Images without Horizon')
clims(2,:) = caxis;


minc = min(clims(:,1));
maxc = max(clims(:,2));

for sb = 3:4
    subplot(2,2,sb)
    caxis([minc maxc])
end

%%

f = fspecial('gaussian',[256,256],24);

horizon_allfixations = imfilter(horizon_allfixations,f);
no_horizon_allfixations = imfilter(no_horizon_allfixations,f);
horizon_BCRW = imfilter(horizon_BCRW,f);
no_horizon_BCRW = imfilter(no_horizon_BCRW,f);

clims = NaN(2,2);

figure
subplot(2,2,1)
imagesc(horizon_allfixations./num_horizon)
colormap('jet')
axis off
axis equal
title('Fixation PDF: Images with Horizon')
clims(1,:) = caxis;

subplot(2,2,2)
imagesc(no_horizon_allfixations./num_no_horizon)
colormap('jet')
axis off
axis equal
title('Fixation PDF: Images without Horizon')
clims(2,:) = caxis;

minc = min(clims(:,1));
maxc = max(clims(:,2));

for sb = 1:2
    subplot(2,2,sb)
    caxis([minc maxc])
end

%%

clims = NaN(2,2);

subplot(2,2,3)
imagesc(horizon_BCRW./num_horizon)
colormap('jet')
axis off
axis equal
title('BCRW: Images with Horizon')
clims(1,:) = caxis;


subplot(2,2,4)
imagesc(no_horizon_BCRW./num_no_horizon)
colormap('jet')
axis off
axis equal
title('BCRW: Images without Horizon')
clims(2,:) = caxis;


minc = min(clims(:,1));
maxc = max(clims(:,2));

for sb = 3:4
    subplot(2,2,sb)
    caxis([minc maxc])
end



%% Determine if there's a seperation in salience and iamge intensity values at fixation locations

horizon_salience_values = [];
horizon_imageI_values = [];
no_horizon_salience_values = [];
no_horizon_imageI_values = [];

for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    
    matfiles = what;
    statfiles = zeros(1,length(tags));
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'FixationStatistics.mat'));
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    statfiles(ii) = i;
                end
            end
        end
    end
    
    for stat = 1:4
        load(matfiles.mat{statfiles(stat)})
        combineddata = shuffunshuffdata{2}{3};
        
        imgI = shuffunshuffdata{2}{3}(:,3,1,1:minlen);
        imgI = reshape(imgI,size(imgI,1),size(imgI,4));
        sal =  shuffunshuffdata{2}{3}(:,1,1,1:minlen);
        sal = reshape(sal,size(sal,1),size(sal,4));
        for i = 1:36;
            if all_horizons{SET}(i) == 1 %has a horizon
                horizon_salience_values = [horizon_salience_values; sal(i,:)];
                horizon_imageI_values = [horizon_imageI_values; imgI(i,:)];
            else
                no_horizon_salience_values = [no_horizon_salience_values; sal(i,:)];
                no_horizon_imageI_values = [no_horizon_imageI_values; imgI(i,:)];
            end
        end
    end
end

figure
subplot(1,2,1)
plot(nanmean(horizon_salience_values))
hold on
plot(nanmean(no_horizon_salience_values),'r')
xlabel('Ordinal Fixation #')
ylabel('Normalized Salience')
legend('Horizon','No Horizon')

subplot(1,2,2)
plot(nanmean(horizon_imageI_values))
hold on
plot(nanmean(no_horizon_imageI_values),'r')
xlabel('Ordinal Fixation #')
ylabel('Normalized ImageI')
legend('Horizon','No Horizon')

%% How is the fit for images with and without Horizon 

horizon_ROC = [];
no_horizon_ROC = [];

horizon_KL= [];
no_horizon_KL = [];
for SET = 1:length(image_sets);
    
    for img = 1:36
        img_ind = 36*(SET-1)+img;
        if all_horizons{SET}(img) == 1 %has a horizon
            horizon_ROC = [horizon_ROC ROC{SET}(:,img)];
            horizon_KL = [horizon_KL; allKL(img_ind,:)];
        else
            no_horizon_ROC = [no_horizon_ROC ROC{SET}(:,img)];
            no_horizon_KL = [no_horizon_KL; allKL(img_ind,:)];
        end
    end
end

p_vals = NaN(1,3);
for i = 1:3
    [~,p_vals(i)] = kstest2(no_horizon_ROC(i,:),horizon_ROC(i,:));
end

figure
b = bar([nanmean(horizon_ROC'); nanmean(no_horizon_ROC')]');
hold on
errorb([nanmean(horizon_ROC'); nanmean(no_horizon_ROC')]',...
    [nanstd(horizon_ROC')./sqrt(size(horizon_ROC,2)) ; nanstd(no_horizon_ROC')./sqrt(size(no_horizon_ROC,2))]')
for c = 1:3
    text(c,0.77,['p = ' num2str(p_vals(c),'%1.1e\n')],'HorizontalAlignment','center')
end
hold off
set(gca,'XTickLabel',{'Fixation vs Salience','Fixation vs BCRW',...
    'Fixation vs Image Intensity'})
legend('Horizon','No Horizon')
yl = ylim;
ylim([0.5 yl(2)])
%%

p_vals = NaN(1,4);
for i = 1:4
    [~,p_vals(i)] = kstest2(no_horizon_KL(:,i),horizon_KL(:,i));
end


figure
b = bar([nanmean(horizon_KL); nanmean(no_horizon_KL)]');
hold on
errorb([nanmean(horizon_KL); nanmean(no_horizon_KL)]',...
    [nanstd(horizon_KL)./sqrt(size(horizon_KL,1)) ; nanstd(no_horizon_KL)./sqrt(size(no_horizon_KL,1))]')
for c = 1:4
    text(c,7.5,['p = ' num2str(p_vals(c),'%1.1e\n')],'HorizontalAlignment','center')
end
hold off
set(gca,'XTickLabel',{'Fixation vs Salience','Fixation vs BCRW',...
    'Fixation vs Image Intensity','Salainece  vs Image Intensity'})

legend('Horizon','No Horizon')
yl = ylim;



