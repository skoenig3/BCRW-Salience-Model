% Get Vivian's Eye Data from Nik.itm
file1 = 'S:\cortex data\vivian\PW121218.3'; %nik1.tim
file2 = 'S:\cortex data\vivian\PW130104.3'; %nik2.itm
file3 = 'S:\cortex data\vivian\PW121220.2'; %nik3.itm
slash = strfind(file1,'\');
slash = slash(end);

[fixationstats] = ClusterFixation(file1);
save([file1(slash+1:end-2) '_' file1(end) '-fixations'],'fixationstats');
[fixationstats] = ClusterFixation(file2);
save([file2(slash+1:end-2) '_' file2(end) '-fixations'],'fixationstats');
[fixationstats] = ClusterFixation(file3);
save([file3(slash+1:end-2) '_' file3(end) '-fixations'],'fixationstats');
%%

itm1 = 'S:\Niklas\Nik1.itm';
itm2 = 'S:\Niklas\Nik2.itm';
itm3 = 'S:\Niklas\Nik3.itm';
images = [];
h = fopen(itm1);
line = fgetl(h);
while line ~= - 1
   line = fgetl(h);
   str = strfind(line,'C:\Nik');
   if ~isempty(str);
       str = line(str+8:end);
       str2 = strfind(str,'.');
      images = [images; str2double(str(1:str2-1))];
   end
end
fclose(h);
[uniqueimages1,novelcnd1,~] = unique(images,'stable');

images = [];
h = fopen(itm2);
line = fgetl(h);
while line ~= - 1
   line = fgetl(h);
   str = strfind(line,'C:\Nik');
   if ~isempty(str);
       str = line(str+8:end);
       str2 = strfind(str,'.');
      images = [images; str2double(str(1:str2-1))];
   end
end
fclose(h);
[uniqueimages2,novelcnd2,~] = unique(images,'stable');

images = [];
h = fopen(itm3);
line = fgetl(h);
while line ~= - 1
   line = fgetl(h);
   str = strfind(line,'C:\Nik');
   if ~isempty(str);
       str = line(str+8:end);
       str2 = strfind(str,'.');
      images = [images; str2double(str(1:str2-1))];
   end
end
fclose(h);
[uniqueimages3,novelcnd3,~] = unique(images,'stable');
%%
imageX = 800;
imageY = 600;
samprate = 5;
PLOTOPTIONS = 'none';

getViewingBehaviorCFNik('PW121218_3-fixations',samprate,imageX,imageY,PLOTOPTIONS,novelcnd1)
getViewingBehaviorCFNik('PW130104_3-fixations',samprate,imageX,imageY,PLOTOPTIONS,novelcnd2)
getViewingBehaviorCFNik('PW121220_2-fixations',samprate,imageX,imageY,PLOTOPTIONS,novelcnd3)
%%
%---[6] Combine Viewing Behavior by Monkey---%
tags = {'PW'};
medianfix = NaN(1,length(tags));
mediansac = NaN(1,length(tags));

matfiles = what;
statfiles = [];
for i = 1:length(matfiles.mat);
    str = strfind(matfiles.mat{i},'ViewingBehavior');
    if ~isempty(str)
        for ii = 1:length(tags);
            strt = strfind(matfiles.mat{i},tags{ii});
            if ~isempty(strt)
                load(matfiles.mat{i},'avgfixprofile','avgsacprofile');
                medianfix(ii) = size(avgfixprofile,2);
                mediansac(ii) = size(avgsacprofile,2);
            end
        end
    end
end
medianfix = round(nanmedian(medianfix));
mediansac = round(nanmedian(mediansac));

allview = cell(1,length(tags));
for i = 1:length(tags)
    allview{i}.densitymap = zeros(600,800);
    allview{i}.allfixations = [];
    allview{i}.allsaccades = [];
    allview{i}.persistence =[];
    allview{i}.anglebtwfix = [];
    allview{i}.sacangle_2fix = [];
    allview{i}.distanceprofile = [];
    allview{i}.distbtwnfix = [];
    allview{i}.fixduration = [];
    allview{i}.sacangle = [];
    allview{i}.sacdist = [];
    allview{i}.sacduration = [];
    allview{i}.timebtwfix = [];
end

matfiles = what;
statfiles = [];
for i = 1:length(matfiles.mat);
    str = strfind(matfiles.mat{i},'ViewingBehavior');
    if ~isempty(str)
        for ii = 1:length(tags);
            strt = strfind(matfiles.mat{i},tags{ii});
            if ~isempty(strt)
                load(matfiles.mat{i});
                if size(allfixations,2) == medianfix(ii);
                    timewarp = 1:size(allfixations,2);
                else
                    timewarp = round(linspace(1,size(allfixations,2),medianfix(ii)));
                end
                distanceprofile.fix = distanceprofile.fix(:,timewarp);
                persistence.fix = persistence.fix(:,timewarp);
                persistence.fix = persistence.fix(:,6:end-5);
                distanceprofile.fix = distanceprofile.fix(:,6:end-5);
                allview{ii}.allfixations = [allview{ii}.allfixations;...
                    allfixations(:,timewarp,:)];
                if size(allsaccades,2) == mediansac(ii);
                    timewarp = 1:size(allsaccades,2);
                else
                    timewarp = round(linspace(1,size(allsaccades,2),mediansac(ii)));
                end
                allview{ii}.allsaccades = [allview{ii}.allsaccades;...
                    allsaccades(:,timewarp,:)];
                distanceprofile.sac = distanceprofile.sac(:,timewarp);
                distanceprofile.sac = distanceprofile.sac(:,6:end-5);
                persistence.sac = persistence.sac(:,timewarp);
                persistence.sac = persistence.sac(:,6:end-5);
                allview{ii}.persistence = [ allview{ii}.persistence;
                    [persistence.sac persistence.fix]];
                allview{ii}.anglebtwfix = [allview{ii}.anglebtwfix;anglebtwfix];
                allview{ii}.sacangle_2fix = [allview{ii}.sacangle_2fix;...
                    sacangle_2fix];
                allview{ii}.densitymap = allview{ii}.densitymap+densitymap;
                allview{ii}.distanceprofile = [allview{ii}.distanceprofile;
                    [distanceprofile.sac distanceprofile.fix]];
                allview{ii}.distbtwnfix = [allview{ii}.distbtwnfix;distbtwnfix];
                allview{ii}.fixduration = [allview{ii}.fixduration;...
                    fixduration];
                allview{ii}.sacangle = [allview{ii}.sacangle;sacangle];
                allview{ii}.sacdist = [allview{ii}.sacdist;sacdist];
                allview{ii}.sacduration = [allview{ii}.sacduration;sacduration];
                allview{ii}.timebtwfix = [allview{ii}.timebtwfix;...
                    timebtwfix];
                allview{ii}.mediansac = mediansac(ii)-10;
                allview{ii}.medianfix = medianfix(ii)-10;
            end
        end
    end
end
clearvars -except allview tags

SAMPRATE = 5;
n = (-180:180)*pi/180;
variables = {'Dist','vel','accel','rot'};
f = fspecial('gaussian',[100,100],50);
graphnum = gcf;
if graphnum == 1
    graphnum = 0;
end
for i = 1:length(tags)
    [allprobanglebtwfix] = hist(allview{i}.anglebtwfix(~isnan(allview{i}.anglebtwfix)),360);
    allprobanglebtwfix = [allprobanglebtwfix(36:-1:1) allprobanglebtwfix allprobanglebtwfix(end:-1:end-36)];
    allprobanglebtwfix = filtfilt(1/6*ones(1,6),1,allprobanglebtwfix);
    allprobanglebtwfix = allprobanglebtwfix(37:end-37);
    allprobanglebtwfix = allprobanglebtwfix/sum(allprobanglebtwfix);
    allprobanglebtwfix = [allprobanglebtwfix allprobanglebtwfix(1)];
    
    figure(graphnum+1)
    title('Distribution of angles between fixations')
    polar(n,allprobanglebtwfix)
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    
    [allprobsacangle] = hist(allview{i}.sacangle(~isnan(allview{i}.sacangle)),360);
    allprobsacangle = [allprobsacangle(36:-1:1) allprobsacangle allprobsacangle(end:-1:end-36)];
    allprobsacangle = filtfilt(1/6*ones(1,6),1,allprobsacangle);
    allprobsacangle = allprobsacangle(37:end-37);
    allprobsacangle = allprobsacangle/sum(allprobsacangle);
    allprobsacangle = [allprobsacangle allprobsacangle(1)];
    
    figure(graphnum+2)
    title('Distribution of angles leaving a fixation')
    polar(n,allprobsacangle)
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    
    %---Stats by fixation---%
    allStatsbyfixation{i}.fixatoinspertrial = reshape(sum(~isnan(allview{i}.fixduration(1:end-1,:)),2),[],66);
    allStatsbyfixation{i}.meanfixationduration = SAMPRATE*nanmean(allview{i}.fixduration);
    allStatsbyfixation{i}.stdfixationduration = SAMPRATE*nanstd(allview{i}.fixduration);
    allStatsbyfixation{i}.numfix = sum(~isnan(allview{i}.fixduration));
    allStatsbyfixation{i}.meansacdistance = nanmean(allview{i}.sacdist);
    allStatsbyfixation{i}.stdsacdistance = nanstd(allview{i}.sacdist);
    allStatsbyfixation{i}.numsacs = sum(~isnan(allview{i}.sacdist));
    
    figure(graphnum+3)
    title('Distribution of Fixation Durations')
    fixduration = allview{i}.fixduration(~isnan(allview{i}.fixduration))*SAMPRATE;
    fixduration(fixduration > 500) = [];
    hist(fixduration,95)
    xlabel('Time (ms)')
    
    figure(graphnum+4)
    title('Distribution of Distances between Fixations')
    hist(allview{i}.distbtwnfix(~isnan(allview{i}.distbtwnfix)),100)
    xlabel('Distance (Pixels)')
    
    figure(graphnum+5)
    title('Fixation Rate')
    hist(1000./(allview{i}.timebtwfix(~isnan(allview{i}.timebtwfix))*SAMPRATE),100)
    xlabel('Hz')
    
    figure(graphnum+6)
    title('Distribution of Saccade Durations')
    sacduration = allview{i}.sacduration(~isnan(allview{i}.sacduration))*SAMPRATE;
    sacduration(sacduration > 150) = [];
    hist(sacduration,25)
    xlabel('Time (ms)')
    
    figure(graphnum+7)
    title('Distribution of Saccade Distances')
    sacdist = allview{i}.sacdist(~isnan(allview{i}.sacdist));
    sacdist(sacdist > 800) = [];
    hist(sacdist,100)
    xlabel('Distance (Pixels)')
    
    figure(graphnum+8)
    title('Correlation between saccade duration and distance')
    plot(allview{i}.sacdist(~isnan(allview{i}.sacdist)),...
        allview{i}.sacduration(~isnan(allview{i}.sacduration))*SAMPRATE,...
        '*','markersize',2)
    xlabel('Distance (pixels)')
    ylabel('Duration (ms)')
    xlim([0 1000])
    ylim([0 200])
    
    figure(graphnum+9)
    title('Probability Distribution of Fixations')
    densitymap = allview{i}.densitymap;
    densitymap = imfilter(densitymap,f);
    densitymap = densitymap./sum(sum(densitymap));
    imagesc(densitymap)
    axis off
    
    figure(graphnum+10)
    title('Fixation Statistics by Trial Number, Fixation Number, or Saccade Number')
    hold all
    errorbar(mean(allStatsbyfixation{i}.fixatoinspertrial),...
        std(allStatsbyfixation{i}.fixatoinspertrial)/...
        size(allStatsbyfixation{i}.fixatoinspertrial,1));
    errorbar(allStatsbyfixation{i}.meanfixationduration,...
        allStatsbyfixation{i}.stdfixationduration./sqrt(allStatsbyfixation{i}.numfix))
    xl = find(allStatsbyfixation{i}.numfix < 50);
    xlim([0 xl(1)])
    ylim([0 400])
    errorbar(allStatsbyfixation{i}.meansacdistance,...
        allStatsbyfixation{i}.stdsacdistance./sqrt(allStatsbyfixation{i}.numsacs))
    if i == 1;
        ylabel('Number of Fixations, Fixation Duration (ms),Saccade Distance (Pixels)')
        set(get(gca,'YLabel'),'Position',[-7 -0.3 0])
    end
    if i == 2
        legend('# Fixations','Fixation Duration','Saccade Distance','Location',...
            'NorthEastOutside')
    end
    if i > 2
        xlabel('Trial Number, Fixation Number, or Saccade Number')
    end
    
    avgfixation= mean(allview{i}.allfixations,1);
    fixlen = size(avgfixation,2);
    avgfixprofile = zeros(size(avgfixation));
    for ii = 1:size(avgfixation,3);
        avgfixprofile(:,:,ii) = filtfilt(1/3*ones(1,3),1,avgfixation(:,:,ii));
        avgfixprofile(:,:,ii) = avgfixprofile(:,:,ii) - min(avgfixprofile(:,:,ii));
        avgfixprofile(:,:,ii) = avgfixprofile(:,:,ii)/max(avgfixprofile(:,:,ii));
    end
    avgsaccade= mean(allview{i}.allsaccades,1);
    saclen = size(avgsaccade,2);
    avgsacprofile = zeros(size(avgsaccade));
    for ii = 1:size(avgsaccade,3);
        avgsacprofile(:,:,ii) = filtfilt(1/3*ones(1,3),1,avgsaccade(:,:,ii));
        avgsacprofile(:,:,ii) =  avgsacprofile(:,:,ii) - min(avgsacprofile(:,:,ii));
        avgsacprofile(:,:,ii) = avgsacprofile(:,:,ii)/max(avgsacprofile(:,:,ii));
    end
    
    figure(graphnum+11)
    title('Average-Smoothed Fixation Profile by Parameter')
    hold all
    h = area(5:fixlen-5,ones(1,fixlen-9));
    set(h,'FaceColor',[.75 .75 .75])
    set(h,'EdgeColor','none')
    for ii =  1:size(avgfixprofile,3);
        plot(avgfixprofile(:,:,ii),'linewidth',2)
    end
    hold off
    xlim([1 fixlen])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
    if i == 1
        ylabel('Normalized Value')
        set(get(gca,'YLabel'),'Position',[-2 -0.2 0])
    end
    if i == 2
        legend([{'fixation'} variables],'Location','NorthEastOutside');
    end
    if i > 2
        xlabel('Warped Time')
    end
    
    figure(graphnum+12)
    title('Average-Smoothed Saccade Profile by Parameter')
    hold all
    h1 = area(1:5,ones(1,5));
    set(h1,'FaceColor',[.75 .75 .75])
    set(h1,'EdgeColor','none')
    h2 = area(saclen-4:saclen,ones(1,5));
    set(h2,'FaceColor',[.75 .75 .75])
    set(h2,'EdgeColor','none')
    for ii = 1:size(avgsacprofile,3)
        p(ii) = plot(avgsacprofile(:,:,ii),'linewidth',2);
    end
    hold off
    xlim([1 saclen])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
    if i == 1
        ylabel('Normalized Value')
        set(get(gca,'YLabel'),'Position',[-2 -0.2 0])
    end
    if i == 2
        legend([h1 p],[{'fixation'} variables],'Location','NorthEastOutside');
    end
    if i > 2
        xlabel('Warped Time')
    end
    
    figure(graphnum+13)
    title('Probability of Saccade Angle Changeing > 45 Degrees')
    hold on
    h  = area(allview{i}.mediansac+1:size(allview{i}.persistence,2),...
        ones(1,size(allview{i}.persistence,2)-allview{i}.mediansac));
    set(h,'FaceColor',[.75 .75 .75])
    set(h,'EdgeColor','none')
    p = plot(mean(allview{i}.persistence));
    hold off
    xlim([1 size(allview{i}.persistence,2)])
    if i == 1
        ylabel('Probability of Saccade Angle Changeing > 45 Degrees')
        set(get(gca,'YLabel'),'Position',[-3 -0.3 0])
    end
    if i == 2
        legend([h p],{'fixation','persistence'},'Location','NorthEastOutside');
    end
    if i > 2 
        xlabel('Warped Time')
    end
    
    figure(graphnum+14)
    plot(nanmean(allview{i}.distanceprofile))
    title(tags{i})
    xlabel('Warped Time')
    ylabel('Distance (pixels)')
    
    [allprobangle2fix] = hist(allview{i}.sacangle_2fix(~isnan(allview{i}.sacangle_2fix)),360);
    allprobangle2fix = [allprobangle2fix(36:-1:1) allprobangle2fix allprobangle2fix(end:-1:end-36)];
    allprobangle2fix = filtfilt(1/6*ones(1,6),1,allprobangle2fix);
    allprobangle2fix = allprobangle2fix(37:end-37);
    allprobangle2fix = allprobangle2fix/sum(allprobangle2fix);
    allprobangle2fix = [allprobangle2fix allprobangle2fix(1)];
    
    figure(graphnum+15)
    title('Distribution of saccade angles entering fixations')
    polar(n,allprobangle2fix)
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
end

screen_size = get(0, 'ScreenSize');
figdir = 'C:\Users\skoenig\Documents\MATLAB\PW Behavioral Data Distributions';
figuretitles = {
    'Distribution of Angles between fixations';
    'Distribution of Angles Leaving a fixation';
    'Distribution of Fixation Durations';
    'Distribution of Distances Between Fixations';
    'Fixation (Saccade) Rate';
    'Distribution of Saccade Durations';
    'Distribution of Saccade Distances';
    'Correlation between Saccade Duration and Saccade Distance';
    'Fixation and Saccade Statistics by Fixation or Saccade Number';
    '2-D Fixation PDF'
    'Average-Smoothed Fixation Profile';
    'Average-Smoothed Saccade Profile';
    'Persistence Profile';
    'Saccade and Fixation Distance Profiles';
    'Distribution of angles entering a fixation';
    };

for ff = 1:15
    figure(graphnum+ff)
    set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    saveas(gcf,[figdir figuretitles{ff}])
    print(gcf,'-r300','-djpeg',[figdir figuretitles{ff}])
end

allviewvariables = {
    'tags: subject names';
    'allview: all combined data by subject';
    'allview.densitymap: positions of fixations';
    'allview.allfixations: fixation profile by parameters warped twice to median length';
    'allview.allsaccades: saccade profile by parameters warped twice to median length';
    'allview.persistence: persistence/probability of eye movement changing >45 angle double warped';
    'allview.anglebtwfix: angles between fixations by fixation';
    'allview.sacangle_2fix: saccade angles entering a fixation';
    'allview.distanceprofile: velocity of eye movements for saccade + subsequent fixation';
    'allview.distbtwnfix: distance between fixations by fixation number';
    'allview.fixduration: fixation duration by fixation number';
    'allview.sacangle: angle of saccade leaving a fixation by saccade number';
    'allview.sacdist: distance of saccade by saccade number';
    'allview.timebtwfix: time between fixations by fixation number';
    };

save('CombinedViewingBehavior-Vivian.mat','tags','allview',...
    'allStatsbyfixation','allviewvariables');