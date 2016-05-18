data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';

scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP'};%,'TT','JN','IW'};

transitionthreshold = 45; %counting as a turn/rotation

medianfix = [46 44 44 44];
mediansac = [21 21 20 20];


vel = NaN(40000,100);
pos_change = NaN(40000,100);
fixdurs = NaN(40000,100);
sacdurs = NaN(40000,100);
fixind = 1;
sacind = 1;
sacind2 = 1;
index = 1;

sac_dir = NaN(40000,100);

sac_persistence = NaN(40000,mediansac(1));
fix_persistence = NaN(40000,medianfix(1));

for imset = 1:length(image_sets);
    dirName = [scm_image_dir image_sets{imset}];
    disp(['Image set-' num2str(image_sets{imset})])
    cd(dirName)
    
    for img = 1:36
        for t = 1:length(tags)
            
            load(['BCRW 2 IOR TAU 17\' tags{t} '-' num2str(img) '-BCRW.mat'],'fixations','fixationstats')
            
            fixlen = medianfix(t); 
            saclen = mediansac(t); 
            for i = 1:20:length(fixationtimes);
                fixdurs(index,1:length(fixationstats{i}.fixationtimes)) = fixationstats{i}.fixationtimes(2,:)-fixationstats{i}.fixationtimes(1,:)+1;
                sacdurs(index,1:length(fixationstats{i}.saccadetimes)) = fixationstats{i}.saccadetimes(2,:)-fixationstats{i}.saccadetimes(1,:)+1;
                index = index+1;
                                
                fixationtimes = fixationstats{i}.fixationtimes;
                saccadetimes = fixationstats{i}.saccadetimes; 
                xy = fixationstats{i}.XY;
                
                for f = 1:size(fixationtimes,2)
                    velx = diff(xy(1,fixationtimes(1,f):fixationtimes(2,f)));
                    vely = diff(xy(2,fixationtimes(1,f):fixationtimes(2,f)));
                    v = sqrt(velx.^2+vely.^2);
                    if length(v) > 100
                        v = v(1:100);
                    end
                    vel(fixind,1:length(v)) = v;
                    
                    
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
                    
                    transitions = transitions(timewarp);
                    
                    if length(fixx) == 1
                        disp('why fix?')
                        continue
                    end
                                        
                    fix_persistence(fixind,:) = transitions;
    
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
                for f = 1:20:size(saccadetimes,2)
                    sacx = xy(1,saccadetimes(1,f):saccadetimes(2,f));
                    sacy = xy(2,saccadetimes(1,f):saccadetimes(2,f));
                    
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
                    
                    if  saclen == length(transitions);
                        timewarp = 1:length(transitions);
                    else
                        timewarp = round(linspace(1,length(transitions),saclen));
                    end
                    
                    transitions = transitions(timewarp);
                    
                    if length(sacx) == 1
                        disp('why sac?')
                        continue
                    end
                    
                    sacangle(f) = atan2(diff(sacy(1:2)),diff(sacx(1:2))); %initial saccade angle
                    
                    sac_persistence(sacind2,:) = transitions;
                    sacind2 = sacind2+1;
                end
                sac_dir(sacind,1:length(sacangle)) = sacangle;
                sacind = sacind +1;
            end
        end
    end
end
fixdur = 5*fixdur(1:end);
sacdur = 5*sacdur(1:end);
%% Histograms for fixation and saccade durations
figure
subplot(1,2,1)
fixdur(fixdur > 500) = [];
hist(fixdur,100);

subplot(1,2,2)
sacdur(sacdur > 100) = [];
hist(sacdur,15);

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
polarplot(n,allprobsacangle)

%% persistence
figure
plot(1-[nanmean(sac_persistence(:,5:end-6)) nanmean(fix_persistence(:,5:end-6))])
hold on
plot([11.5 11.5],[0 1],'k--')
hold off
