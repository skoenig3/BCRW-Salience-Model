scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;
binsize = 25; %24 pixels per dva but using 25 cuz divides nicely into 600 & 800
f = fspecial('gaussian',[155,155],24); %1 dva smoothing
distance_threshold = 48;

IOR_dir = 'BCRW IOR TAU Simulations Saccadic Momentum\';
IOR_taus = [0 35 30 25 20 17 14 12 9 7 5 3 1];


%define diagonals for 1-back,2-back,3-back
id0 =  eye(40);
id0 =  [id0(2:end,:); zeros(1,40)];
id1 = eye(40);
id1 = [id1(3:end,:); zeros(2,40)];
id2 = eye(40);
id2 = [id2(4:end,:); zeros(3,40)];
id3 = eye(40);
id3 = [id3(5:end,:); zeros(4,40)];
id4 = eye(40);
id4 = [id4(6:end,:); zeros(5,40)];


all_recurrence = cell(1,length(IOR_taus));
all_corms = cell(1,length(IOR_taus));
all_rates = cell(1,length(IOR_taus));
for it = 1:length(IOR_taus)
    all_recurrence{it} =zeros(40,40);
    all_corms{it} = NaN(1,288);
    all_rates{it} = NaN(1,288); 
end

%%
for imset = 1:length(image_sets);
    disp(['Image set-' num2str(image_sets{imset})])
    dirName = [scm_image_dir image_sets{imset} '\'];
    
    matfiles = what(dirName);
    
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
    
    for IOR = 1:length(IOR_taus);
        for img = 1:36;
            load([dirName num2str(img) '-saliencemap.mat'],'fullmap');
            allfixations = zeros(imageY,imageX);
            allBCRW = zeros(imageY,imageX);
            BCRW_at_fixations = [];
            BCRW_at_random = [];
            index = (imset-1)*36+img;
            
            sal_at_fixations = NaN(400,50);
            count = 1;
            for t = 1:length(tags)
                if IOR_taus(IOR) == 0
                    load([dirName IOR_dir tags{t} '-' num2str(img) '-' num2str(IOR_taus(IOR)) '-BCRW-SM.mat'])
                else
                    load([dirName IOR_dir tags{t} '-' num2str(img) '-' num2str(1/IOR_taus(IOR)) '-BCRW-SM.mat'])
                end
                 BCRWfixations = fixationtimes;
                
                %---Calculate Salience at BCRW Fixations---%
                corm = NaN(1,size(BCRWfixations,1),1);
                rec_rate = NaN(1,size(BCRWfixations,1),1);
                for ii = 1:size(BCRWfixations,1);
                    recurence_map = zeros(40,40);
                    tind = find(BCRWfixations(ii,:,1) > 0);
                    if length(tind) > 50;
                        tind = tind(1:50);
                    end
                    fixxy =  [BCRWfixations(ii,tind,1); BCRWfixations(ii,tind,2)];
                    if size(fixxy,2) > 40
                        fixxy = fixxy(:,1:40);
                    end
                   
                    [rate,recurence_map,this_corm] = ...
                        calculate_auto_recurrence(fixxy,40,distance_threshold);
                    
%                     remove 0-backs
%                     rm = recurence_map;
%                     rm(id0 == 1) = 0;
%                     rm(id0' == 1) = 0;
%                     rm(tril(rm) == 1) = 0;
%                     if sum(sum(rm)) > 3;
%                         [xcm,ycm] = centroid(rm);
%                         corm(count) = xcm-ycm;
%                         count = count+1;
%                     end
                    
                    if ~isnan(this_corm)
                        corm(count) = this_corm;
                        rec_rate(count) = rate;
                        count = count+1;
                    end
                    
                    
                    all_recurrence{IOR} = all_recurrence{IOR}+recurence_map;
                end
            end
            all_corms{IOR}(index) = nanmean(corm);
            all_rates{IOR}(index) = nanmean(rec_rate);
        end
    end
end
%%
normalized_all_recurrence = all_recurrence;
for it = 1:length(IOR_taus)
    rm =  all_recurrence{it};
    for r = 1:length(rm);
        fix_count = rm(r,r);
        for i = 1:r-1
            rm(r,i) = rm(r,i)/fix_count;
            rm(i,r) = rm(i,r)/fix_count;
        end
        rm(r,r) = 1;
    end
    rm(find(eye(size(rm)))) = NaN;
%     rm(id0 == 1) = NaN;
%     rm(id0' == 1) = NaN;
%     rm(id1 == 1) = NaN;
%     rm(id1' == 1) = NaN;
%     rm(id2 == 1) = NaN;
%     rm(id2' == 1) = NaN;
%     rm(id3 == 1) = NaN;
%     rm(id3' == 1) = NaN;
%     rm(id4 == 1) = NaN;
%     rm(id4' == 1) = NaN;
    normalized_all_recurrence{it} = rm;
end

figure
for it = 2:2:length(IOR_taus)
    subplot(2,3,(it+0)/2)
    h = imagesc(100*normalized_all_recurrence{it});
    axis xy, set(h,'alphadata',~isnan(normalized_all_recurrence{it}));
    colorbar
    box off
    axis square
    xlabel('Ordinal Fixation #')
    ylabel('Ordinal Fixation #')
    xlim([0.5 35])
    ylim([0.5 35])
    title(['\tau_{IOR} = ' num2str(IOR_taus(it))])
    colormap('viridis')
end

%%
all_coms = [];
all_its = [];
all_coms_means = [];
all_coms_stds = [];
xlabl = {};
for IOR = 1:length(IOR_taus)
    all_coms = [all_coms  all_corms{IOR}(1,:)];%(1,:)-all_corms{IOR}(2,:)];
    all_its = [all_its IOR*ones(1,size( all_corms{IOR},2))];
    xlabl = [xlabl {['\tau 1/' num2str(IOR_taus(IOR))]}];
    
    all_coms_means = [all_coms_means mean(all_corms{IOR}(1,:))];%-all_corms{IOR}(2,:))];
    all_coms_stds = [all_coms_stds std(all_corms{IOR}(1,:))];%-all_corms{IOR}(2,:))];
end

[P,T,STATS,TERMS] = anovan(all_coms',all_its');
COMPARISON = multcompare(STATS);
pk = kruskalwallis(all_coms,all_its');

figure
hold on
bar(all_coms_means)
errorb(all_coms_means,all_coms_stds./sqrt(size(all_corms{1},2)))
hold off
set(gca,'XTick',[1:length(IOR_taus)])
set(gca,'XTickLabel',xlabl)
ylabel('CORM,Local->Global')
ylim([3 4])

%%
all_rats = [];
all_its = [];
all_rats_means = [];
all_rats_stds = [];
xlabl = {};
for IOR = 1:length(IOR_taus)
    all_rats = [all_rats  all_rates{IOR}(1,:)];%(1,:)-all_rates{IOR}(2,:)];
    all_its = [all_its IOR*ones(1,size( all_rates{IOR},2))];
    xlabl = [xlabl {['\tau 1/' num2str(IOR_taus(IOR))]}];
    
    all_rats_means = [all_rats_means mean(all_rates{IOR}(1,:))];%-all_rates{IOR}(2,:))];
    all_rats_stds = [all_rats_stds std(all_rates{IOR}(1,:))];%-all_rates{IOR}(2,:))];
end

[P,T,STATS,TERMS] = anovan(all_rats',all_its');
COMPARISON = multcompare(STATS);
pk = kruskalwallis(all_rats,all_its');
%%
figure
bar(all_rats_means)
hold on
errorb(all_rats_means,all_rats_stds./sqrt(size(all_rates{1},2)))
hold off
%%
set(gca,'XTick',[1:length(IOR_taus)])
set(gca,'XTickLabel',xlabl)
ylabel('Recurrence Rate (%)')
ylim([3 4])