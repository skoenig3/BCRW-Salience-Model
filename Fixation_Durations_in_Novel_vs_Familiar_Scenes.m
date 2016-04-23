scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009'};%,...
%     'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','JN','IW','TT'};
fixdur = NaN(288,75,4,2);
numfixes = NaN(288,4,2);
imageX = 800;
imageY = 600;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            for ii = 1:length(tags);
                strt = strfind(matfiles.mat{i},tags{ii});
                if ~isempty(strt)
                    load(matfiles.mat{i})
                    for cndlop = 1:2:length(fixationstats)
                        fixationtimes = fixationstats{cndlop}.fixationtimes;
                        if ~isempty(fixationtimes)
                            fixations = fixationstats{cndlop}.fixations;
                            if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                                    fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                                fixationtimes(:,1) = [];
                            end
                            numfixes((SET-1)*36+(cndlop+1)/2,ii,1) = size(fixationtimes,2);
                            fixdur((SET-1)*36+(cndlop+1)/2,1:size(fixationtimes,2),ii,1)  = ...
                                diff(fixationtimes,1)'+1;
                        end
                    end
                    for cndlop = 2:2:length(fixationstats)
                        if strcmp(trialtype(cndlop),'f')
                            fixationtimes = fixationstats{cndlop}.fixationtimes;
                            fixations = fixationstats{cndlop}.fixations;
                            if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                                    fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                                fixationtimes(:,1) = [];
                            end
                            numfixes((SET-1)*36+cndlop/2,ii,2) = size(fixationtimes,2);
                            fixdur((SET-1)*36+cndlop/2,1:size(fixationtimes,2),ii,2)  = ...
                                diff(fixationtimes,1)'+1;
                        end
                    end
                end
            end
        end
    end
end
%% Novel vs Familiar for JN MP IW
novelnumfixes = numfixes(:,1:3,1);
novelnumfixes = novelnumfixes(1:end);
famnumfixes = numfixes(:,1:3,2);
famnumfixes = famnumfixes(1:end);

[h,p] = ttest2(novelnumfixes,famnumfixes)
%% Novel vs Familiar for TT
novelnumfixes = numfixes(:,4,1);
novelnumfixes = novelnumfixes(1:end);
famnumfixes = numfixes(:,4,2);
famnumfixes = famnumfixes(1:end);

[h,p] = ttest2(novelnumfixes,famnumfixes)
%%
nov = fixdur(:,1:35,1:3,1);
nov = 5*[nov(:,:,1);nov(:,:,2);nov(:,:,3)];
rep = fixdur(:,1:35,1:3,2);
rep = 5*[rep(:,:,1);rep(:,:,2);rep(:,:,3)];
nov(nov > 660) = NaN;
rep(rep > 660) = NaN;
figure
hold on
errorbar(nanmean(nov),nanstd(nov)./sqrt(sum(~isnan(nov))))
errorbar(nanmean(rep),nanstd(rep)./sqrt(sum(~isnan(rep))),'r')
hold off
legend('Novel','Familiar')
title('Fixation Duration for MP, JN, IW, and TT pre-lesion')
xlabel('Fixation Number')
ylabel('Fixation Duration (ms)')
%%
novdur = fixdur(:,1:35,1:4,1);
novdur = novdur(1:end); 
famdur = fixdur(:,1:35,1:4,2);
famdur = famdur(1:end);
[h, p] = kstest2(novdur,famdur)
%%
figure
plot(5*median(nanmedian(fixdur(:,1:35,4,1)),3))
hold on
plot(5*median(nanmedian(fixdur(1:144,1:35,4,2)),3),'r')
plot(5*median(nanmedian(fixdur(145:end,1:35,4,2)),3),'g')
hold off
legend('Novel','Familiar pre-lesion','Familiar post-lesion')
title('Fixation Duration for TT')
xlabel('Fixation Number')
ylabel('Fixation Duration (ms)')
%%
novdur = fixdur(:,1:35,4,1);
novdur = novdur(1:end); 
famdur = fixdur(:,1:35,4,2);
famdur = famdur(1:end);
[h, p] = kstest2(novdur,famdur)
