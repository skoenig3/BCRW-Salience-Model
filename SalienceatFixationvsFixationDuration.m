% Salience at Fixation Location vs Fixation Duration
scm_image_dir = 'C:\Users\GOD-OF-ChAOS\Documents\MATLAB\Buffalo Lab-Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
PLOTOPTIONS = 'none';
imageX = 800; imageY = 600;
saldur = NaN(45000,2);
fixcount = 1;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    disp(SETNUM)
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
        load(matfiles.mat{eyefile})
        
        matfiles = what;
        saliencemapfiles = [NaN;NaN];
        for i = 1:length(matfiles.mat);
            str = strfind(matfiles.mat{i},'saliencemap');
            if ~isempty(str)
                dash = strfind(matfiles.mat{i},'-');
                saliencemapfiles = [saliencemapfiles [i;str2num(matfiles.mat{i}(1:dash(1)-1))]];
            end
        end
        saliencemapfiles(:,1) = [];
        [~,si] = sort(saliencemapfiles(2,:));
        saliencemapfiles = si;
        
        for cndlop=1:2:length(fixationstats)
            reindexed = (cndlop+1)/2;
            load(matfiles.mat{saliencemapfiles(reindexed)})
            
            fixationtimes = fixationstats{cndlop}.fixationtimes;
            fixations = fixationstats{cndlop}.fixations;
            xy = fixationstats{cndlop}.XY;
            if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                    fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                fixations(:,1) = [];
                fixationtimes(:,1) = [];
            end
           
            for i = 1:size(fixations,2);
                spot = ceil(fixations(:,i));
                spot(2) = imageY-spot(2);
                spot(spot < 1) = 1;
                spot(1,spot(1) > imageX) = imageX;
                spot(2,spot(2) > imageY) = imageY;
                
                saldur(fixcount,1) = fullmap(spot(2),spot(1));
                saldur(fixcount,2) = fixationtimes(2,i)-fixationtimes(1,i);
                fixcount = fixcount+1;
            end
        end
    end
end
%%
nan = find(isnan(saldur(:,1)));
saldur(nan,:) = [];
toolong = find(saldur(:,2) > 100);
saldur(toolong,:) = [];