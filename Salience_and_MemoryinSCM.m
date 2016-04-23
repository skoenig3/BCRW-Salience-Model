% scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
% image_sets = {'Set006','Set007','Set008','Set009'};
% screen_size = get(0, 'ScreenSize');
% imageX = 800; imageY = 600;
%
% load('C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\CriticalRegions.mat');
% salience = cell(2,2);
% numfixations = cell(2,2);
% fixationnumber = cell(2,2);
% for SET = 1:length(image_sets);
%     SETNUM = image_sets{SET};
%     cd([scm_image_dir SETNUM])
%     matfiles = what;
%     eyedatafiles = [];
%     for i = 1:length(matfiles.mat);
%         str = strfind(matfiles.mat{i},'fixation');
%         if ~isempty(str)
%             eyedatafiles = [eyedatafiles i];
%         end
%     end
%     for eyefile = eyedatafiles;
%         load(matfiles.mat{eyefile});
%         for i = 2:2:length(fixationstats);
%             novfixations = fixationstats{i-1}.fixations;
%             if ~isempty(novfixations);
%                 novfixations(:,2) = imageY-novfixations(:,2);
%                 novfixations(novfixations < 1) = 1;
%                 novfixations(novfixations(:,1) > imageX,1) = imageX;
%                 novfixations(novfixations(:,2) > imageY,2) = imageY;
%                 repfixations = fixationstats{i}.fixations;
%                 repfixations(:,2) = imageY-repfixations(:,2);
%                 repfixations(repfixations < 1) = 1;
%                 repfixations(repfixations(:,1) > imageX,1) = imageX;
%                 repfixations(repfixations(:,2) > imageY,2) = imageY;
%                 if strcmp(trialtype(i),'r')
%                     ROI = allROIs{SET,1}(i/2,:);
%                     ROIfix = find(...
%                         novfixations(1,:) > ROI(1) & novfixations(1,:) < ROI(2) & ...
%                         novfixations(2,:) > ROI(3) & novfixations(2,:) < ROI(4));
%                     if ~isempty(ROIfix)
%                         load([num2str(i/2) '-saliencemap.mat'],'fullmap');
%                         sal = 0;
%                         for f = 1:length(ROIfix);
%                             sal = sal+fullmap(round(novfixations(2,ROIfix(f))),...
%                                 round(novfixations(1,ROIfix(f))));
%                         end
%                         sal = sal/f;
%                         fixationnumber{1,1} = [fixationnumber{1,1} ROIfix(1)];
%                         salience{1,1} = [salience{1,1} sal];
%                         numfixations{1,1} = [numfixations{1,1} f];
%                     else
%                         fixationnumber{1,1} = [fixationnumber{1,1} NaN];
%                         salience{1,1} = [salience{1,1} NaN];
%                         numfixations{1,1} = [numfixations{1,1} 0];
%                     end
%                     ROIfix = find(...
%                         repfixations(1,:) > ROI(1) & repfixations(1,:) < ROI(2) & ...
%                         repfixations(2,:) > ROI(3) & repfixations(2,:) < ROI(4));
%                     if ~isempty(ROIfix)
%                         load([num2str(i/2) '-saliencemap.mat'],'fullmap');
%                         sal = 0;
%                         for f = 1:length(ROIfix);
%                             sal = sal+fullmap(round(repfixations(2,ROIfix(f))),...
%                                 round(repfixations(1,ROIfix(f))));
%                         end
%                         sal = sal/f;
%                         salience{1,2} = [salience{1,2} sal];
%                         fixationnumber{1,2} = [fixationnumber{1,2} ROIfix(1)];
%                         numfixations{1,2} = [numfixations{1,2} f];
%                     else
%                         salience{1,2} = [salience{1,2} NaN];
%                         fixationnumber{1,2} = [fixationnumber{1,2} NaN];
%                         numfixations{1,2} = [numfixations{1,2} 0];
%                     end
%                 elseif strcmp(trialtype(i),'m')
%                     ROI = allROIs{SET,3}(i/2,:);
%                     ROIfix = find(...
%                         novfixations(1,:) > ROI(1) & novfixations(1,:) < ROI(2) & ...
%                         novfixations(2,:) > ROI(3) & novfixations(2,:) < ROI(4));
%                     if ~isempty(ROIfix)
%                         load([num2str(i/2) '-saliencemap.mat'],'fullmap');
%                         sal = 0;
%                         for f = 1:length(ROIfix);
%                             sal = sal+fullmap(round(novfixations(2,ROIfix(f))),...
%                                 round(novfixations(1,ROIfix(f))));
%                         end
%                         sal = sal/f;
%                         fixationnumber{2,1} = [fixationnumber{2,1} ROIfix(1)];
%                         salience{2,1} = [salience{2,1} sal];
%                         numfixations{2,1} = [numfixations{2,1} f];
%                     else
%                         fixationnumber{2,1} = [fixationnumber{2,1} NaN];
%                         salience{2,1} = [salience{2,1} NaN];
%                         numfixations{2,1} = [numfixations{2,1} NaNrre];
%                     end
%                     ROI = allROIs{SET,2}(i/2,:);
%                     ROIfix = find(...
%                         repfixations(1,:) > ROI(1) & repfixations(1,:) < ROI(2) & ...
%                         repfixations(2,:) > ROI(3) & repfixations(2,:) < ROI(4));
%                     if ~isempty(ROIfix)
%                         load([num2str(i/2) '-saliencemap.mat'],'fullmap');
%                         sal = 0;
%                         for f = 1:length(ROIfix);
%                             sal = sal+fullmap(round(repfixations(2,ROIfix(f))),...
%                                 round(repfixations(1,ROIfix(f))));
%                         end
%                         sal = sal/f;
%                         salience{2,2} = [salience{2,2} sal];
%                         fixationnumber{2,2} = [fixationnumber{2,2} ROIfix(1)];
%                         numfixations{2,2} = [numfixations{2,2} f];
%                     else
%                         salience{2,2} = [salience{2,2} NaN];
%                         fixationnumber{2,2} = [fixationnumber{2,2} NaN];
%                         numfixations{2,2} = [numfixations{2,2} NaN];
%                     end
%                 end
%             end
%         end
%     end
% end
%%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009'};
screen_size = get(0, 'ScreenSize');
imageX = 800; imageY = 600;

load('C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\CriticalRegions.mat');
salience = cell(2,2);
time = cell(2,2);
firsttime = cell(2,2);
area = cell(2,2);
imgtype = 'om';
missed = cell(3,2,2);
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    for eyefile = eyedatafiles;
        load(matfiles.mat{eyefile});
        for i = 2:2:length(fixationstats);
            novfixations = fixationstats{i-1}.fixations;
            novtimes =  5*fixationstats{i-1}.fixationtimes;
            novfixations(2,:) = imageY-novfixations(2,:);
            repfixations = fixationstats{i}.fixations;
            reptimes =  5*fixationstats{i}.fixationtimes;
            repfixations(2,:) = imageY-repfixations(2,:);
            if strcmp(trialtype(i),'r') || strcmp(trialtype(i),'m')
                salmap1 = [num2str(i/2) '-saliencemap.mat'];
                if strcmp(trialtype(i),'r')
                    salmap2 = [num2str(i/2) 'o-saliencemap.mat'];
                    ROI1 = allROIs{SET,1}(i/2,:);
                    ROI1(1) = ROI1(1)-24;
                    ROI1(2) = ROI1(2)+24;
                    ROI1(3) = ROI1(3)-24;
                    ROI1(4) = ROI1(4)+24;
                    ROI1(ROI1 < 1) = 1;
                    ROI1(ROI1 > imageX) = imageX;
                    if ROI1(3) > imageY; ROI1(3) = imageY;end
                    if ROI1(4) > imageY; ROI1(4) = imageY;end
                    ROI2 =ROI1;
                    type = 1;
                elseif strcmp(trialtype(i),'m');
                    salmap2 = [num2str(i/2) 'm-saliencemap.mat'];
                    ROI1 = allROIs{SET,3}(i/2,:);
                    ROI1(1) = ROI1(1)-24;
                    ROI1(2) = ROI1(2)+24;
                    ROI1(3) = ROI1(3)-24;
                    ROI1(4) = ROI1(4)+24;
                    ROI1(ROI1 < 1) = 1;
                    ROI1(ROI1 > imageX) = imageX;
                    if ROI1(3) > imageY; ROI1(3) = imageY;end
                    if ROI1(4) > imageY; ROI1(4) = imageY;end
                    ROI2 = allROIs{SET,2}(i/2,:);
                    ROI2(1) = ROI2(1)-24;
                    ROI2(2) = ROI2(2)+24;
                    ROI2(3) = ROI2(3)-24;
                    ROI2(4) = ROI2(4)+24;
                    ROI2(ROI2 < 1) = 1;
                    ROI2(ROI2 > imageX) = imageX;
                    if ROI2(3) > imageY; ROI2(3) = imageY;end
                    if ROI2(4) > imageY; ROI2(4) = imageY;end
                    type = 2;
                end
                
                
                
                %                         figure
                %                         set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
                %                         subplot(1,2,1)
                %                         imshow(imread([num2str(i/2) '.bmp']))
                %                         hold on
                %                         plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
                %                             [ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
                %                         plot(novfixations(1,:),novfixations(2,:))
                %                         hold off
                %                         axis off
                %                         subplot(1,2,2)
                %                         imshow(imread([num2str(i/2) imgtype(type) '.bmp']));
                %                         hold on
                %                         plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
                %                             [ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
                %                         plot(repfixations(1,:),repfixations(2,:))
                %                         hold off
                %                         axis off
                
                
                ROIfix1 = find(...
                    novfixations(1,:) > ROI1(1) & novfixations(1,:) < ROI1(2) & ...
                    novfixations(2,:) > ROI1(3) & novfixations(2,:) < ROI1(4));
                ROIfix2 = find(...
                    repfixations(1,:) > ROI2(1) & repfixations(1,:) < ROI2(2) & ...
                    repfixations(2,:) > ROI2(3) & repfixations(2,:) < ROI2(4));
                
                %                         subplot(1,2,1)
                %                         hold on
                %                         for r = 1:length(ROIfix1);
                %                             plot(novfixations(1,ROIfix1(r)),novfixations(2,ROIfix1(r)),'*m');
                %                         end
                %                         hold off
                %                         subplot(1,2,2)
                %                         hold on
                %                         for r = 1:length(ROIfix2);
                %                             plot(repfixations(1,ROIfix2(r)),repfixations(2,ROIfix2(r)),'*m');
                %                         end
                %                         hold off
                %
                
                if ~isempty(ROIfix1) &&  ~isempty(ROIfix2);
                    load(salmap1,'fullmap');
                    sal = fullmap(ROI1(3):ROI1(4),ROI1(1):ROI1(2));
                    sal(sal < mean2(sal) + 1/2*std2(sal)) = NaN;
                    salience{1,type} = [salience{1,type} nanmedian(nanmedian(sal))];
                    time{1,type} = [time{1,type} length(ROIfix1)/size(novfixations,2)];
                    firsttime{1,type} = [firsttime{1,type} novtimes(1,ROIfix1(1))];
                    area{1,type} = [area{1,type} (ROI1(1)-ROI1(2))*(ROI1(3)-ROI1(4))];
                    load(salmap2,'fullmap')
                    sal = fullmap(ROI2(3):ROI2(4),ROI2(1):ROI2(2));
                    sal(sal < mean2(sal) + 1/2*std2(sal)) = NaN;
                    salience{2,type} = [salience{2,type} nanmedian(nanmedian(sal))];
                    time{2,type} = [time{2,type} length(ROIfix2)/size(repfixations,2)];
                    firsttime{2,type} = [firsttime{2,type} reptimes(1,ROIfix2(1))];
                    area{2,type} = [area{2,type} (ROI2(1)-ROI2(2))*(ROI2(3)-ROI2(4))];
                else
                    load(salmap1,'fullmap');
                    sal = fullmap(ROI1(3):ROI1(4),ROI1(1):ROI1(2));
                    sal(sal < mean2(sal) + 1/2*std2(sal)) = NaN;
                    missed{1,type,1} = [missed{1,type,1} nanmedian(nanmedian(sal))];
                    if ~isempty(ROIfix1)
                        missed{2,type,1} = [missed{2,type,1} length(ROIfix1)/size(novfixations,2)];
                        missed{3,type,1}= [missed{3,type,1} novtimes(1,ROIfix1(1))];
                    else
                        missed{2,type,1} = [missed{2,type,1} NaN];
                        missed{3,type,1} = [missed{3,type,1} NaN];
                        
                    end
                    
                    load(salmap2,'fullmap');
                    sal = fullmap(ROI2(3):ROI2(4),ROI2(1):ROI2(2));
                    sal(sal < mean2(sal) + 1/2*std2(sal)) = NaN;
                    missed{1,type,2} = [missed{1,type,2}  nanmedian(nanmedian(sal))];
                    if ~isempty(ROIfix2)
                        missed{2,type,2} =[missed{2,type,2} length(ROIfix2)/size(repfixations,2)];
                        missed{3,type,2}= [missed{3,type,2} reptimes(1,ROIfix2(1))];
                    else
                        missed{2,type,2} = [missed{2,type,2} NaN];
                        missed{3,type,2} = [missed{3,type,2} NaN];
                    end
                end
                %                         close
            end
        end
    end
end

toolong = find(firsttime{2,1} > 2500);
firsttime{2,1}(toolong) = [];
time{2,1}(toolong) = [];
salience{2,1}(toolong) = [];
firsttime{1,1}(toolong) = [];
time{1,1}(toolong) = [];
salience{1,1}(toolong) = [];
area{1,1}(toolong) = [];
area{2,1}(toolong) = [];
toolong = find(firsttime{2,2} > 2500);
firsttime{2,2}(toolong) = [];
time{2,2}(toolong) = [];
salience{2,2}(toolong) = [];
firsttime{1,2}(toolong) = [];
time{1,2}(toolong) = [];
salience{1,2}(toolong) = [];
area{1,2}(toolong) = [];
area{2,2}(toolong) = [];
%% For Missed objects before manipulation-is it a difference in salience?
novsalo = missed{1,1,1}; novsalo(isnan(missed{3,1,1})) = [];
novsalm = missed{1,2,1}; novsalm(isnan(missed{3,2,1})) = [];
[~,psalo] = kstest2(salience{1,1},novsalo)
[~,psalm] = kstest2(salience{1,2},novsalm)
%% For Missed objects with manipulation-is it a difference in salience?
novsalo = missed{1,1,2}; novsalo(isnan(missed{3,1,2})) = [];
novsalm = missed{1,2,2}; novsalm(isnan(missed{3,2,2})) = [];

[~,psalo] = kstest2(salience{2,1},novsalo)
[~,psalm] = kstest2(salience{2,2},novsalm)
%% For objects that were missed in novel presentation, what was the salience values when found during manipulation
spontfound = missed{1,1,2}(isnan(missed{2,1,1}));
[~,pspontfoundo] = kstest2(salience{2,1},spontfound)

spontfound = missed{1,2,2}(isnan(missed{2,2,1}));
[~,pspontfoundm] = kstest2(salience{2,2},spontfound)
%% For objects that were missed in novel presentation, what was their previous salience value before manipulation
spontfound = missed{1,1,1}(isnan(missed{2,1,1}));
[~,pspontfoundo] = kstest2(salience{1,1},spontfound)

spontfound = missed{1,2,1}(isnan(missed{2,2,1}));
[~,pspontfoundm] = kstest2(salience{1,2},spontfound)
%% Perecent of Fixations in ROI
[~,pnumfixeso] = kstest2(time{1,1},time{2,1})
[~,pnumfixesm] = kstest2(time{1,2},time{2,2})

means = [mean(time{1,1}) mean(time{2,1});
    mean(time{1,2}) mean(time{2,2})];
%% Perecent of Fixations in ROI normalized by area
means = [mean(time{1,1}./area{1,1}) mean(time{2,1}./area{2,1});
    mean(time{1,2}./area{1,2}) mean(time{2,2}./area{2,2})];

[~,pnumfixeso] = kstest2(time{1,1}./area{1,1},time{2,1}./area{2,1})
[~,pnumfixesm] = kstest2(time{1,2}./area{1,2},time{2,2}./area{2,2})
%% Time to first fixation [onset] in ROI
[~,pfix1o] = kstest2(firsttime{1,1},firsttime{2,1})
[~,pfix1m] = kstest2(firsttime{1,2},firsttime{2,2})

means = [mean(firsttime{1,1}) mean(firsttime{2,1});
    mean(firsttime{1,2}) mean(firsttime{2,2})];

%% Time to first fixation [onset] in ROI normalized by area
% may not need to normalize by area because not percentage of fixations
% though area may make some changes more aparent if larger
means = [mean(firsttime{1,1}./area{1,1}) mean(firsttime{2,1}./area{2,1});
    mean(firsttime{1,2}./area{1,2}) mean(firsttime{2,2}./area{2,2})];

[~,pfix1o] = kstest2(firsttime{1,1}./area{1,1},firsttime{2,1}./area{2,1})
[~,pfix1m] = kstest2(firsttime{1,2}./area{1,2},firsttime{2,2}./area{2,2})
%%
