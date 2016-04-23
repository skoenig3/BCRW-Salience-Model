% % Written by Seth Koenig 9/3/15
% % Code losely based on PrePostLesionAnalysis.m & PrePostLesionAnalysisANOVA.m,
% % but largely rewritten to streamline.
% % Code determines significant sources of various in behavioral measures
% % such as salience at fixation locations and fixation durations.
% % Importanlty the code can be used to determine if there's differences
% % between viewing behavior across SCM vs SCME sets and pre-post lesion for
% % TT.
%
% %---ANOVA for Salience and Image intesnity @ Fixation Locations---%
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';

scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};

load([data_dir 'CombinedSalienceStatistics.mat'],'alldata') %salience at fixation locations
salience_at_fixation = reshape(alldata(:,1,1,:),[size(alldata,1) size(alldata,4)]); %
salience_at_fixation = salience_at_fixation(:,1:35); %organized every 4th row per monkey, and median num fix is 35
imgintensity_at_fixation = reshape(alldata(:,3,1,:),[size(alldata,1) size(alldata,4)]); %image intensity at fixation locations
imgintensity_at_fixation = imgintensity_at_fixation(:,1:35); %organized every 4th row per monkey, and median num fix is 35


names = {};
sets = [];
salvals = [];
ivals = [];
group = [];
fixationnumber = [];
for t = 1:length(tags) %for analysis of TT pre-post lesion set t = 2 only
    for SET = 1:length(image_sets)
        data = salience_at_fixation(t:4:end,:); %salience for that monkey for all sets
        data = data(36*(SET-1)+1:SET*36,:); %data just for that set

        idata = imgintensity_at_fixation(t:4:end,:); %salience for that monkey for all sets
        idata = idata(36*(SET-1)+1:SET*36,:); %data just for that set
        for ii = 1:size(data,1);
            data2 = data(ii,:);
            data2(isnan(data2)) = [];

            idata2 = idata(ii,:);
            idata2(isnan(idata2)) = [];
            if  ~isempty(data2)
                %monkey name
                nms = repmat(tags{t},[length(data2),1]);
                nms = cellstr(nms);
                names = [names;nms];
                %set number
                sts = repmat(SET,[length(data2),1]);
                sets = [sets;sts];
                %set_group
                if SET <= 4
                    group = [group;ones(length(data2),1)];
                else
                    group = [group;2*ones(length(data2),1)];
                end
                %salience values
                salvals = [salvals;data2'];
                %image intensity values
                ivals = [ivals;idata2'];
                %fixation number
                fixationnumber = [fixationnumber; [1:length(data2)]'];
            end
        end
    end
end
psal = anovan(salvals,{names,sets,fixationnumber});
psal2= anovan(salvals,{group}); %run seperately since sets are not independent of group (aka df = 0)
psalI = anovan(ivals,{names,sets,fixationnumber});
psalI2= anovan(ivals,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%
%---ANOVA for Viewing behavior---%
%not all variable have the same number of data points e.g # of saccades ~= # of fixations all of the time,
%and certain varialbe have 1 less data point/trial because they calculate the difference
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';

scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;

load([data_dir 'CombinedViewingBehavior.mat'],'allview') %all viewing behavior statistics

%---Analysis of Saccade Stuff---%
%Stat variables
sacangle_2fix = []; %saccade angle entering fixation
sacangle = []; %saccade angle leaving a fixation
sacdist = []; %saccade arc length
sacamplitude = []; %saccae amplitude
sacduration = []; %saccade duration

% Factor variables
names = []; %monkey name
sets = []; %image set
group = []; %SCM vs SCME
imgnum = [];%image so ignoring sets
fixationnumber = [];%fixation number or saccade number
for t = 1:length(tags) %for analysis of TT pre-post lesion set t = 2 only
    for SET = 1:length(image_sets)
        for ii = (SET-1)*36+1:SET*36;
            
            %---Stats stuff---%
            %saccade angle leaving fixation
            data = allview{t}.sacangle(ii,:);
            data(isnan(data)) = [];
            sacangle = [sacangle; data'];
            
            %saccade angle entering fixation
            data = allview{t}.sacangle_2fix(ii,:);
            data(isnan(data)) = [];
            sacangle_2fix = [sacangle_2fix; data'];
            
            %saccade arclength
            data = allview{t}.sacdist(ii,:);
            data(isnan(data)) = [];
            sacdist = [sacdist; data'];
            
            %saccade amplitude
            data = allview{t}.sacamplitude(ii,:);
            data(isnan(data)) = [];
            sacamplitude = [sacamplitude; data'];
            
            %saccade duration
            data = allview{t}.sacduration(ii,:);
            data(isnan(data)) = [];
            sacduration = [sacduration; data'];
            
            %--Factor stuff---%
            %monkey name
            nms = t*ones(length(data),1);
            names = [names;nms];
            %set number
            sts = SET*ones(length(data),1);
            sets = [sets;sts];
            %set_group
            if SET <= 4
                group = [group;ones(length(data),1)];
            else
                group = [group;2*ones(length(data),1)];
            end
            %fixation number
            fixationnumber = [fixationnumber; [1:length(data)]'];
            %image
            imgnum = [imgnum; ii*ones(length(data),1)];
        end
    end
end

%%
psacangle = anovan(sacangle,{names,fixationnumber,imgnum});
psacangle2 = anovan(sacangle,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%
psacangle_2fix = anovan(sacangle_2fix,{names,sets,fixationnumber,imgnum});
psacangle_2fix2 = anovan(sacangle_2fix,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%
psacdist = anovan(sacdist,{names,sets,fixationnumber,imgnum});
psacdist2 = anovan(sacdist,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%
psacamplitude = anovan(sacamplitude,{names,sets,fixationnumber,imgnum});
psacamplitude2 = anovan(sacamplitude,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%
psacduration = anovan(sacduration,{names,sets,fixationnumber,imgnum});
psacduration2 = anovan(sacduration,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%

%---Analysis of Fixation Stuff---%
%number of fixations don't necessarily equal number of saccades
%and for some reason time betwen fix is different probably because
%doesn't calculate first fixation the same way (see next section)
%Stat variables
distbtwnfix = [];%distance between fixations
fixduration = []; %fixation durations
anglebtwfix = []; %angle between fixations

% Factor variables
names = []; %monkey name
sets = []; %image set
group = []; %SCM vs SCME
imgnum = [];%image so ignoring sets
fixationnumber = [];%fixation number or saccade number
for t = 1:length(tags) %for analysis of TT pre-post lesion set t = 2 only
    for SET = 1:length(image_sets)
        for ii = (SET-1)*36+1:SET*36;
        
            
            %angle between fixations
            data = allview{t}.anglebtwfix(ii,:);
            data(isnan(data)) = [];
            anglebtwfix = [anglebtwfix; data'];
            
            %distance between fixations
            data = allview{t}.distbtwnfix(ii,:);
            data(isnan(data)) = [];
            distbtwnfix = [distbtwnfix; data'];
            
            %fixation duration
            data = allview{t}.fixduration(ii,:);
            data(isnan(data)) = [];
            fixduration = [fixduration; data'];
            
            %--Factor stuff---%
            %monkey name
            nms = t*ones(length(data),1);
            names = [names;nms];
            %set number
            sts = SET*ones(length(data),1);
            sets = [sets;sts];
            %set_group
            if SET <= 4
                group = [group;ones(length(data),1)];
            else
                group = [group;2*ones(length(data),1)];
            end
            %fixation number
            fixationnumber = [fixationnumber; [1:length(data)]'];
            %image
            imgnum = [imgnum; ii*ones(length(data),1)];
        end
    end
end

%%
pfixduration = anovan(fixduration,{names,sets,fixationnumber,imgnum});
pfixduration2 = anovan(fixduration,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%
pdistbtwnfix = anovan(distbtwnfix,{names,sets,fixationnumber,imgnum});
pdistbtwnfix2 = anovan(distbtwnfix,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%
panglebtwfix = anovan(anglebtwfix,{names,sets,fixationnumber,imgnum});
panglebtwfix2 = anovan(anglebtwfix,{group}); %run seperately since sets are not independent of group (aka df = 0)
%%

%Stat variables
timebtwfix = []; %time between fixations or 1/saccade rate

% Factor variables
names = []; %monkey name
sets = []; %image set
group = []; %SCM vs SCME
imgnum = [];%image so ignoring sets
fixationnumber = [];%fixation number or saccade number
for t = 1:length(tags) %for analysis of TT pre-post lesion set t = 2 only
    for SET = 1:length(image_sets)
        for ii = (SET-1)*36+1:SET*36;
        
            %---Stats stuff---%      
            %time between fixations i.e. 1/saccade rate
            data = allview{t}.timebtwfix(ii,:);
            data(isnan(data)) = [];
            timebtwfix = [timebtwfix; data'];
            
            %--Factor stuff---%
            %monkey name
            nms = t*ones(length(data),1);
            names = [names;nms];
            %set number
            sts = SET*ones(length(data),1);
            sets = [sets;sts];
            %set_group
            if SET <= 4
                group = [group;ones(length(data),1)];
            else
                group = [group;2*ones(length(data),1)];
            end
            %fixation number
            fixationnumber = [fixationnumber; [1:length(data)]'];
            %image
            imgnum = [imgnum; ii*ones(length(data),1)];
        end
    end
end

ptimebtwfix = anovan(timebtwfix,{names,sets,fixationnumber,imgnum});
ptimebtwfix2 = anovan(timebtwfix,{group}); %run seperately since sets are not independent of group (aka df = 0)
    
