scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
imageX = 800;
imageY = 600;



% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];

entropyvalues = NaN(1,288);
salience_entropyvalues = NaN(1,288);
edgevalues =NaN(1,288);
img_count = 1;

for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    
    disp(['Set #' num2str(SET)])
    for img = 1:36
        
        
        imgr = imread([num2str(img) '.bmp']);
        load([num2str(img) '-saliencemap.mat'],'fullmap')
        
        entropyvalues(img_count) = entropy(imgr);%pixel intesnity entropy
        salience_entropyvalues(img_count) = entropy(fullmap);
        xedges = imfilter(imgr,sobelx);
        yedges = imfilter(imgr,sobely);
        edgevalues(img_count) = mean2(xedges+yedges); %edgineess
        
        img_count = img_count+1;
    end
end
