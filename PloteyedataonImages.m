img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
figure_dir1 = 'C:\Users\seth.koenig\Documents\MATLAB\SCM eye traces by monkey\';
figure_dir2 = 'C:\Users\seth.koenig\Documents\MATLAB\SCM eye traces by image\';
for set = 1:4
    cd([img_dir 'SetE00' num2str(set)]);
    list = ls;
    eyedat = [];
    for l = 1:size(list,1);
        if ~isempty(strfind(list(l,:),'fixation.mat'))
            eyedat = [eyedat l];
        end
    end
    
    load(list(eyedat(1),:));
    fix1 = fixationstats;
    trialtype1 = trialtype;
    if length(trialtype1) == 71 %happened with 1 of Irwins files
       disp('Missing file')
       trialtype1(72) ='f';  
       fix1(end+1) = fix1(end);
    end
    
    load(list(eyedat(2),:));
    fix2 = fixationstats;
    trialtype2 = trialtype;
    
    load(list(eyedat(3),:));
    fix3 = fixationstats;
    trialtype3 = trialtype;
    
    try
        load(list(eyedat(4),:));
        fix4 = fixationstats;
        trialtype4 = trialtype;
    catch
        fix4 = [];
    end
    
%     for trial = 1:2:72
%         figure
%         subplot(1,2,1)
%         hold on
%         imshow(imread([num2str((trial+1)/2) '.bmp']))
%         plot(fix1{trial}.XY(1,:),600-fix1{trial}.XY(2,:))
%         hold off
%         title('Novel image')
%         box off
%         
%         subplot(1,2,2)
%         hold on
%         if trialtype1(trial+1) == 'f'
%             imshow(imread([num2str((trial+1)/2) '.bmp']))
%             title('Familiar')
%         elseif trialtype1(trial+1) == 'r'
%             imshow(imread([num2str((trial+1)/2) 'o.bmp']))
%             title('Replaced')
%         elseif trialtype1(trial+1) == 'm'
%             imshow(imread([num2str((trial+1)/2) 'm.bmp']))
%             title('Moved')
%         end
%         if ~isempty(fix1{trial+1})
%             plot(fix1{trial+1}.XY(1,:),600-fix1{trial+1}.XY(2,:))
%         end
%         hold off
%         box off
%         
%         subtitle(['Monkey: ' list(eyedat(1),1:2) '  Set ' num2str(set) 'Image ' num2str((trial+1)/2)])
%         save_and_close_fig(figure_dir1,[list(eyedat(1),1:2) '_Set' num2str(set) '_img' num2str((trial+1)/2)])
%     end
%     
%     for trial = 1:2:72
%         figure
%         subplot(1,2,1)
%         hold on
%         imshow(imread([num2str((trial+1)/2) '.bmp']))
%         plot(fix2{trial}.XY(1,:),600-fix2{trial}.XY(2,:))
%         hold off
%         title('Novel image')
%         box off
%         
%         subplot(1,2,2)
%         hold on
%         if trialtype1(trial+1) == 'f'
%             imshow(imread([num2str((trial+1)/2) '.bmp']))
%             title('Familiar')
%         elseif trialtype1(trial+1) == 'r'
%             imshow(imread([num2str((trial+1)/2) 'o.bmp']))
%             title('Replaced')
%         elseif trialtype1(trial+1) == 'm'
%             imshow(imread([num2str((trial+1)/2) 'm.bmp']))
%             title('Moved')
%         end
%         plot(fix2{trial+1}.XY(1,:),600-fix2{trial+1}.XY(2,:))
%         hold off
%         box off
%         
%         subtitle(['Monkey: ' list(eyedat(2),1:2) '  Set ' num2str(set) ' Image ' num2str((trial+1)/2)])
%         save_and_close_fig(figure_dir1,[list(eyedat(2),1:2) '_Set' num2str(set) '_img' num2str((trial+1)/2)])
%     end
%     
%     for trial = 1:2:72
%         figure
%         subplot(1,2,1)
%         hold on
%         imshow(imread([num2str((trial+1)/2) '.bmp']))
%         plot(fix3{trial}.XY(1,:),600-fix3{trial}.XY(2,:))
%         hold off
%         title('Novel image')
%         box off
%         
%         subplot(1,2,2)
%         hold on
%         if trialtype1(trial+1) == 'f'
%             imshow(imread([num2str((trial+1)/2) '.bmp']))
%             title('Familiar')
%         elseif trialtype1(trial+1) == 'r'
%             imshow(imread([num2str((trial+1)/2) 'o.bmp']))
%             title('Replaced')
%         elseif trialtype1(trial+1) == 'm'
%             imshow(imread([num2str((trial+1)/2) 'm.bmp']))
%             title('Moved')
%         end
%         plot(fix3{trial+1}.XY(1,:),600-fix3{trial+1}.XY(2,:))
%         hold off
%         box off
%         
%         subtitle(['Monkey: ' list(eyedat(3),1:2) '  Set ' num2str(set) 'Image ' num2str((trial+1)/2)])
%         save_and_close_fig(figure_dir1,[list(eyedat(3),1:2) '_Set' num2str(set) '_img' num2str((trial+1)/2)])
%     end
%     
%     for trial = 1:2:72
%         figure
%         subplot(1,2,1)
%         hold on
%         imshow(imread([num2str((trial+1)/2) '.bmp']))
%         plot(fix4{trial}.XY(1,:),600-fix4{trial}.XY(2,:))
%         hold off
%         title('Novel image')
%         box off
%         
%         subplot(1,2,2)
%         hold on
%         if trialtype1(trial+1) == 'f'
%             imshow(imread([num2str((trial+1)/2) '.bmp']))
%             title('Familiar')
%         elseif trialtype1(trial+1) == 'r'
%             imshow(imread([num2str((trial+1)/2) 'o.bmp']))
%             title('Replaced')
%         elseif trialtype1(trial+1) == 'm'
%             imshow(imread([num2str((trial+1)/2) 'm.bmp']))
%             title('Moved')
%         end
%         if ~isempty(fix4{trial+1})
%             plot(fix4{trial+1}.XY(1,:),600-fix4{trial+1}.XY(2,:))
%         end
%         hold off
%         box off
%         
%         subtitle(['Monkey: ' list(eyedat(4),1:2) '  Set ' num2str(set) ' Image ' num2str((trial+1)/2)])
%         save_and_close_fig(figure_dir1,[list(eyedat(4),1:2) '_Set' num2str(set) '_img' num2str((trial+1)/2)])
%     end
    
    for trial = 1:2:72
        figure
        
        subplot(2,2,1)
        hold on
        imshow(imread([num2str((trial+1)/2) '.bmp']))
        plot(fix1{trial}.XY(1,:),600-fix1{trial}.XY(2,:))
        hold off
        title('Novel image')
        box off
        
        subplot(2,2,2)
        hold on
        imshow(imread([num2str((trial+1)/2) '.bmp']))
        plot(fix2{trial}.XY(1,:),600-fix2{trial}.XY(2,:))
        hold off
        title('Novel image')
        box off
        
        subplot(2,2,3)
        hold on
        imshow(imread([num2str((trial+1)/2) '.bmp']))
        plot(fix3{trial}.XY(1,:),600-fix3{trial}.XY(2,:))
        hold off
        title('Novel image')
        box off
        
        subplot(2,2,4)
        hold on
        imshow(imread([num2str((trial+1)/2) '.bmp']))
        plot(fix4{trial}.XY(1,:),600-fix4{trial}.XY(2,:))
        hold off
        title('Novel image')
        box off
        
        subtitle(['Set ' num2str(set) ' Image ' num2str((trial+1)/2)])
        save_and_close_fig(figure_dir2,['Set' num2str(set) '_img' num2str((trial+1)/2)])
    end
end
