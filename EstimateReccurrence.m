% Calculate Recurence

scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
imageX = 800;
imageY = 600;

tags = {'MP','TT','JN','IW'};
recurence_measure = [];
all_recurence_map = zeros(50,50);
count = 1;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'fixation'))
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    eyedatafiles(ii) = i;
                end
            end
        end
    end
    
    for eyefile = eyedatafiles;
        load(matfiles.mat{eyefile})
        for cndlop=1:2:length(fixationstats)
            fixations = fixationstats{cndlop}.fixations;
            if ~isempty(fixations)
                if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                        fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                    fixations(:,1) = [];
                end
                
                recurence_map = zeros(50,50);
                N=size(fixations,2);
                N(N > 50) = 50;
                [x,y]=meshgrid(1:N);
                i=find(ones(N)); %forms pairs except for self-pairing
                i=[x(i), y(i)];
                dist =sqrt((fixations(1,i(:,1))-fixations(1,i(:,2))).^2 +...
                    (fixations(2,i(:,1))-fixations(2,i(:,2))).^2);
                dind = find(dist <= 48);
                for d = 1:length(dind);
                    recurence_map(i(dind(d),1),i(dind(d),2)) = recurence_map(i(dind(d),1),i(dind(d),2))+1;
                end
                recurence_measure(count) = 100*2*sum(sum(triu(recurence_map)))/N/(N-1);
                count = count+1;
                %                 figure
                %                 hold on
                %                 for d = 1:length(dind);
                %                     plot(i(dind(d),1),i(dind(d),2),'.')
                %                 end
                %                 hold off
                %                 xlabel('Ordinal Fixation Number')
                %                 ylabel('Ordinal Fixation Number')
                %                 close
            end
            all_recurence_map = all_recurence_map + recurence_map;
        end
    end
end
%%
figure
imagesc(all_recurence_map(end:-1:1,:))
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')

rm = all_recurence_map;
rm(find(eye(size(rm)))) = 0;
figure
imagesc(rm(end:-1:1,:))
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
%% remove central diagonals
id = eye(size(rm));
id = [id(2:end,:); zeros(1,50)];
rm(find(id)) = 0;

id = eye(size(rm));
id = [zeros(1,50);id(1:end-1,:)];
rm(find(id)) = 0;
imagesc(rm(end:-1:1,:))
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
%%

id = eye(size(rm));
id = [id(3:end,:); zeros(2,50)];
rm(find(id)) = 0;

id = eye(size(rm));
id = [zeros(2,50);id(1:end-2,:)];
rm(find(id)) = 0;
imagesc(rm(end:-1:1,:))
xlabel('Ordinal Fixation Number')
ylabel('Ordinal Fixation Number')
xlim([0.5 35.5]) %only care about the 1st 35 fixations
ylim([15.5 50.5])
axis square
%%
figure
hist(recurence_measure,100);
xlim([0 60])
xlabel('Recurrence in %')
ylabel('Count')
title(['Distribution of Recurrence Measures; mean=' num2str(mean(recurence_measure))])