scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
screen_size = get(0, 'ScreenSize');
imageX = 800; imageY = 600;
% allROIs = cell(8,3);
for SET = 5%:length(image_sets);
    dirName = [scm_image_dir image_sets{SET}];
    cd(dirName)
    
    for i = 1:36;
        for ii = 1:3;
            allROIs{SET,ii} = NaN(36,4);
        end
    end
    for i = 4%1:36
        
        
        img = num2str(i);
        novimg = imread([img '.bmp']);
        replimg = imread([img 'o.bmp']);
        movedimg = imread([img 'm.bmp']);
        
        replacedchange = abs(sum(novimg-replimg,3));
        movedchange = abs(sum(novimg-movedimg,3));
        
        movedchange(movedchange < 5) = 0;
        replacedchange(replacedchange < 5) = 0;
        
        [r, c] = find(replacedchange ~=0);
        
        replacedROI = [min(c) max(c) min(r) max(r)];
        white = [replacedROI(1)-50 replacedROI(2)+50 replacedROI(3)-50 replacedROI(4)+50];
        white(white < 1) = 1;
        if white(2) > imageX; white(2) = imageX; end
        if white(4) > imageY; white(4) = imageY; end
        movedchange(white(3):white(4),white(1):white(2)) = 0;
        [r, c] = find(movedchange ~= 0);
        movedROI = [min(c) max(c) min(r) max(r)];
        
        if ~isempty(movedROI)
            movedchange = abs(sum(novimg-movedimg,3));
            movedchange(movedchange < 5) = 0;
            white = [movedROI(1)-50 movedROI(2)+50 movedROI(3)-50 movedROI(4)+50];
            white(white < 1) = 1;
            if white(2) > imageX; white(2) = imageX; end
            if white(4) > imageY; white(4) = imageY; end
            movedchange(white(3):white(4),white(1):white(2)) = 0;
            [r, c] = find(movedchange ~= 0);
            if ~isempty(r)
                b4movedROI = [min(c) max(c) min(r) max(r)];
            end
            
            allROIs{SET,1}(i,:) = replacedROI;
            allROIs{SET,2}(i,:) = movedROI;
            allROIs{SET,3}(i,:) = b4movedROI;
        end
        
        figure
        subtitle(['Set ' num2str(SET) ' Image # ' num2str(i)])
        set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        subplot(2,1,1)
        imshow(novimg);
        hold on
        if ~isempty(b4movedROI)
            plot([b4movedROI(1) b4movedROI(2) b4movedROI(2) b4movedROI(1) b4movedROI(1)],...
                [b4movedROI(3) b4movedROI(3) b4movedROI(4) b4movedROI(4) b4movedROI(3)],'g');
        end
        hold off
        title('original image')
        
        subplot(2,2,3)
        imshow(replimg);
        title('replaced image')
        hold on
        plot([replacedROI(1) replacedROI(2) replacedROI(2) replacedROI(1) replacedROI(1)],...
            [replacedROI(3) replacedROI(3) replacedROI(4) replacedROI(4) replacedROI(3)],'r');
        hold off
        
        subplot(2,2,4)
        imshow(movedimg);
        title('moved image')
        hold on
        if ~isempty(movedROI)
            plot([movedROI(1) movedROI(2) movedROI(2) movedROI(1) movedROI(1)],...
                [movedROI(3) movedROI(3) movedROI(4) movedROI(4) movedROI(3)],'r');
            plot([b4movedROI(1) b4movedROI(2) b4movedROI(2) b4movedROI(1) b4movedROI(1)],...
                [b4movedROI(3) b4movedROI(3) b4movedROI(4) b4movedROI(4) b4movedROI(3)],'g');
        end
        hold off
        pause
        close
    end
end