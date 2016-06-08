data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';

scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP'};%,'TT','JN','IW'};

transitionthreshold = 45; %counting as a turn/rotation

medianfix = [46 44 44 44];
mediansac = [21 21 20 20];


fixvel = NaN(400000,medianfix(1));
sacvel = NaN(400000,mediansac(1));
pos_change = NaN(400000,100);
fixdurs = NaN(40000,100);
sacdurs = NaN(40000,100);
fixind = 1;
sacind = 1;
sacind2 = 1;
index = 1;

sac_dir = NaN(400000,100);

sac_persistence = NaN(400000,mediansac(1)-10);
fix_persistence = NaN(400000,medianfix(1)-10);

for imset = 1:length(image_sets);
    dirName = [scm_image_dir image_sets{imset}];
    disp(['Image set-' num2str(image_sets{imset})])
    cd(dirName)
    
    for img = 1:36
        for t = 1:length(tags)
            
            load(['BCRW IOR TAU 17\' tags{t} '-' num2str(img) '-BCRW.mat'],'fixations','fixationstats')
            
            fixlen = medianfix(t);
            saclen = mediansac(t);
            for i = 1:10:length(fixationstats);
                fixdurs(index,1:length(fixationstats{i}.fixationtimes)) = fixationstats{i}.fixationtimes(2,:)-fixationstats{i}.fixationtimes(1,:)+1;
                sacdurs(index,1:length(fixationstats{i}.saccadetimes)) = fixationstats{i}.saccadetimes(2,:)-fixationstats{i}.saccadetimes(1,:)+1;
                
                if any(fixdurs(index,:) < 5)
                    disp('what')
                end
                
                if any(sacdurs(index,:) < 2)
                    disp('what')
                end
                index = index+1;
                
                fixationtimes = fixationstats{i}.fixationtimes;
                saccadetimes = fixationstats{i}.saccadetimes;
                xy = fixationstats{i}.XY;
                
                for f = 1:size(fixationtimes,2)
                    velx = diff(xy(1,fixationtimes(1,f):fixationtimes(2,f)));
                    vely = diff(xy(2,fixationtimes(1,f):fixationtimes(2,f)));
                    
                    
                    
                    fixx = xy(1,fixationtimes(1,f):fixationtimes(2,f));
                    fixy = xy(2,fixationtimes(1,f):fixationtimes(2,f));
                    
                    if f == 1
                        x = xy(1,fixationtimes(1,f):fixationtimes(2,f)+7); %7 not 5 because going to index +2
                        y = xy(2,fixationtimes(1,f):fixationtimes(2,f)+7); %7 not 5 because going to index +2
                    elseif f == size(fixationtimes,2)
                        x = xy(1,fixationtimes(1,f)-5:fixationtimes(2,f));
                        y = xy(2,fixationtimes(1,f)-5:fixationtimes(2,f));
                    else
                        x = xy(1,fixationtimes(1,f)-5:fixationtimes(2,f)+7); %7 not 5 because going to index +2
                        y = xy(2,fixationtimes(1,f)-5:fixationtimes(2,f)+7); %7 not 5 because going to index +2
                    end
                    velx = diff(x);
                    vely = diff(y);
                    angle = 180*atan2(vely,velx)/pi;
                    transitions = abs(diff(angle)) > transitionthreshold;
                    dist = zeros(1,length(x)-2);
                    if  fixlen == length(transitions);
                        timewarp = 1:length(transitions);
                    else
                        timewarp = round(linspace(1,length(transitions),fixlen));
                    end
                    
                    v = sqrt(velx.^2+vely.^2);
                    v = v(timewarp);
                    fixvel(fixind,1:length(v)) = v;
                    
                    transitions = transitions(timewarp);
                    
                    if length(fixx) == 1
                        error('why fix?')
                    end
                    
                    fix_persistence(fixind,:) = transitions(5:end-6);
                    
                    x1 = xy(1,fixationtimes(1,f));
                    y1 = xy(2,fixationtimes(1,f));
                    x = xy(1,fixationtimes(1,f):fixationtimes(2,f));
                    y = xy(2,fixationtimes(1,f):fixationtimes(2,f));
                    
                    dist = sqrt((x-x1).^2+(y-y1).^2);
                    
                    if length(dist) > 100;
                        dist = dist(1:100);
                    end
                    pos_change(fixind,1:length(dist)) = dist;
                    
                    fixind = fixind+1;
                end
                
                sacangle = NaN(1,size(saccadetimes,2));
                for f = 1:size(saccadetimes,2)
                    sacx = xy(1,saccadetimes(1,f):saccadetimes(2,f));
                    sacy = xy(2,saccadetimes(1,f):saccadetimes(2,f));
                    
                    velx = diff(sacx);
                    vely = diff(sacy);
                    
                    
                    
                    
                    
                    if f == 1
                        x = xy(1,saccadetimes(1,f):saccadetimes(2,f)+7); %7 not 5 because going to index +2
                        y = xy(2,saccadetimes(1,f):saccadetimes(2,f)+7); %7 not 5 because going to index +2
                    elseif f == size(saccadetimes,2)
                        x = xy(1,saccadetimes(1,f)-5:saccadetimes(2,f));
                        y = xy(2,saccadetimes(1,f)-5:saccadetimes(2,f));
                    else
                        x = xy(1,saccadetimes(1,f)-5:saccadetimes(2,f)+7); %7 not 5 because going to index +2
                        y = xy(2,saccadetimes(1,f)-5:saccadetimes(2,f)+7); %7 not 5 because going to index +2
                    end
                    velx = diff(x);
                    vely = diff(y);
                    angle = 180*atan2(vely,velx)/pi;
                    transitions = abs(diff(angle)) > transitionthreshold;
                    dist = zeros(1,length(x)-2);
                    
                    vel = sqrt(velx.^2+vely.^2);
                    
                    
                    if  saclen == length(transitions);
                        timewarp = 1:length(transitions);
                    else
                        timewarp = round(linspace(1,length(transitions),saclen));
                    end
                    
                    transitions = transitions(timewarp);
                    vel = vel(timewarp);
                    sacvel(sacind2,1:length(vel)) = vel;
                    
                    if length(sacx) == 1
                        error('why sac?')
                    end
                    
                    if f == 1
                        sacx = xy(1,saccadetimes(1,f):saccadetimes(2,f));
                        sacy = xy(2,saccadetimes(1,f):saccadetimes(2,f));
                    else
                        sacx = xy(1,saccadetimes(1,f)-1:saccadetimes(2,f));
                        sacy = xy(2,saccadetimes(1,f)-1:saccadetimes(2,f));
                    end
                    
                    
                    sacangle(f) = atan2(diff(sacy(1:2)),diff(sacx(1:2))); %initial saccade angle
                    
                    sac_persistence(sacind2,:) = transitions(5:end-6);
                    sacind2 = sacind2+1;
                end
                sac_dir(sacind,1:length(sacangle)) = sacangle;
                sacind = sacind +1;
            end
        end
    end
end
fixdurs = 5*fixdurs(1:end);
sacdurs = 5*sacdurs(1:end);
%%
load([data_dir 'CombinedViewingBehavior.mat'],'allview')
sacduration = allview{1}.sacduration;
sacduration = 5*sacduration(1:end);
sacduration(sacduration > 100) = [];
sac = hist(sacduration,10:5:100);

fixduration = allview{1}.fixduration;
fixduration = 5*fixduration(1:end);
fixduration(fixduration > 500) = [];
fix = hist(fixduration,20:10:500);




%% Histograms for fixation and saccade durations
figure

fixdurs(fixdurs > 500) = [];
fx = hist(fixdurs,20:10:500);
subplot(1,2,1)
plot(20:10:500,fx/sum(fx))
hold on
plot(20:10:500,fix/sum(fix),'r')
hold off
ylabel('Probability')
xlabel('Fixation Duration (ms)')
legend('Simulated','Observed')



sacdurs(sacdurs > 100) = [];
sc = hist(sacdurs,10:5:100);
subplot(1,2,2)
plot(10:5:100,sc/sum(sc));
hold on
plot(10:5:100,sac/sum(sac),'r')
hold off
ylabel('Probability')
xlabel('Saccade Duration (ms)')
legend('Simulated','Observed')



%% Polar plot of saccades leaving fixation
sac_dir = sac_dir(1:end);
sac_dir(isnan(sac_dir)) = [];
n = (-180:180)*pi/180;

[allprobsacangle] = hist(sac_dir,360);
allprobsacangle = [allprobsacangle(36:-1:1) allprobsacangle allprobsacangle(end:-1:end-36)];
allprobsacangle = filtfilt(1/6*ones(1,6),1,allprobsacangle);
allprobsacangle = allprobsacangle(37:end-37);
allprobsacangle = allprobsacangle/sum(allprobsacangle);
allprobsacangle = [allprobsacangle allprobsacangle(1)];

figure
polarplot(n,allprobsacangle/sum(allprobsacangle))
hold on

sacangle = allview{1}.sacangle;
sacangle = sacangle(1:end);
sacangle(isnan(sacangle)) = [];

[allprobsacangle] = hist(sacangle,360);
allprobsacangle = [allprobsacangle(36:-1:1) allprobsacangle allprobsacangle(end:-1:end-36)];
allprobsacangle = filtfilt(1/6*ones(1,6),1,allprobsacangle);
allprobsacangle = allprobsacangle(37:end-37);
allprobsacangle = allprobsacangle/sum(allprobsacangle);
allprobsacangle = [allprobsacangle allprobsacangle(1)];

polarplot(n,allprobsacangle/sum(allprobsacangle),'r')
legend('Simulated','Observed')
%% movement velocity
figure
plot([nanmean(sacvel(:,5:end-6)) nanmean(fixvel(:,5:end-6))]/24)
hold on
plot(nanmean(allview{1}.distanceprofile)/24,'r')

xlabel('Warped Time')
ylabel('Distance (dva)')
box off
legend('Simulated','Observed')
%% persistence
figure
plot(1-[nanmean(sac_persistence) nanmean(fix_persistence)])
hold on
plot(1-nanmean(allview{1}.persistence),'r')
plot([11.5 11.5],[0 1],'k--')
hold off
xlim([1 47])
xlabel('Warped Time')
ylabel('Persistence')
legend('Simulated','Observed')

%% Position change for observed data
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP'};%,'TT','JN','IW'};
PLOTOPTIONS = 'none';
imageX = 800;
imageY = 600;
pos_change = NaN(400000,100);
SAMPRATE = 5;
fixind = 1;
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
    
    for i = 1:length(eyedatafiles)
        load(matfiles.mat{eyedatafiles(i)},'fixationstats')
        
        for i = 1:36
            fixationtimes = fixationstats{2*i-1}.fixationtimes;
            xy = fixationstats{2*i-1}.XY;
            
            for f = 1:size(fixationtimes,2)
                x1 = xy(1,fixationtimes(1,f));
                y1 = xy(2,fixationtimes(1,f));
                x = xy(1,fixationtimes(1,f):fixationtimes(2,f));
                y = xy(2,fixationtimes(1,f):fixationtimes(2,f));
                
                dist = sqrt((x-x1).^2+(y-y1).^2);
                
                if length(dist) > 100;
                    dist = dist(1:100);
                end
                pos_change(fixind,1:length(dist)) = dist;
                fixind = fixind+1;
            end
        end
    end
end