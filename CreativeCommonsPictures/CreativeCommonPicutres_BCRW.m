%code written to show visual example of goodness of fit of BCRW. Had to
%download new images with known copyrights specifically public domain or
%creative commons licenses. Eye tracking data collected from JN who was
%originally in the study.
%
% Code parallels RunWholeSalienceModel code

%---[1] Get the salience maps for multiple task sets---%
% image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\CreativeCommonsPictures\Renamed\';
% cd(image_dir)
% dirData = dir(image_dir);
% dirIndex = [dirData.isdir];
% fileList = {dirData(~dirIndex).name}';
% imageindexes = [];
% for i = 1:length(fileList)
%     bmps = strfind(fileList{i},'bmp');
%     if ~isempty(bmps)
%         imageindexes = [imageindexes i];
%     end
% end
% for i = 1:length(imageindexes)
%     imagefile = fileList{imageindexes(i)};
%     getSalienceMap(imagefile)
% end

%%
%---[8] Run BCRW ---%
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\CreativeCommonsPictures\Renamed\';
tags = {'MP','TT','JN','IW'};
Combinedbehaviorfile = ['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets'...
    '\CombinedViewingBehavior.mat'];
load(Combinedbehaviorfile,'allview')
imageX = 800; imageY = 600;
plotoptions.runs = 'none'; %all/none
plotoptions.probdens = 'none';
plotoptions.type = 'sal'; %sal/image
IOR_taus = 1/17;


cd(image_dir)
% 
% for i = 1:54
%     for t = 1:length(tags)
%         disp(['Running ' tags{t} ' on image #' num2str(i)])
%         run_BCRWCF_saveXY(allview{t},[num2str(i) '-saliencemap'],tags{t},imageX,imageY,plotoptions,IOR_taus)
%     end
% end

%%
tags = {'MP','TT','JN','IW'};
f = fspecial('gaussian',[256,256],24);
for i = 1:54
    
    figure
    subplot(1,3,1)
    imshow(imread([num2str(i) '.bmp']))
    axis off
    
    load([num2str(i) '-saliencemap.mat'],'fullmap')
    subplot(1,3,2)
    imagesc(fullmap)
    axis off
    axis equal
    colormap('jet')
    
    allBCRW = zeros(imageY,imageX);
    for t = 1:length(tags)
        load(['BCRW IOR TAU 17\' tags{t} '-' num2str(i) '-BCRW.mat'],'fixations')
        allBCRW = allBCRW+fixations;
    end
    allBCRW = imfilter(allBCRW,f);
    
    subplot(1,3,3)
    imagesc(allBCRW)
    axis off
    axis equal
    colormap('jet')
    
    save_and_close_fig('Figures\',[num2str(i) '_maps'])
end
    