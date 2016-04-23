%% Analyzing Perecent coverage
% how much did the monkeys and BCRW cover an image with their attention.
% Assumes that with 2 dva (radius) of a fixation location is attended to
% since this is the optimal size for the area affected by IOR in the BCRW,
% the approximate distance between prior and return fixations, and roughly
% the size of the fovea.
%% Percent Coverage for Observed data for all time points (even after 10s)
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = ['MP';'TT';'JN';'IW'];
imageX = 800;
imageY = 600;
samprate = 5;

IOR_area = 48;
[rr,cc] = meshgrid(1:imageX,1:imageY);

count = ones(1,4);
percent_coverage = NaN(4,288);
for imset = 1:length(image_sets);
    cd([image_dir '\' image_sets{imset} '\'])
    a = what;
    for aa = 1:size(a.mat,1);
        if ~isempty(strfind(a.mat{aa},'fixation.mat'))
            disp(['Image set-' image_sets{imset} ' File ' a.mat{aa}])
            datafile = a.mat{aa}(1:10);
            datafile(9) = '.';
            for noreason = 1 %just so I can shorten code visually
                %sub function convert raw eye tracking data into x & y coordinates
                [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(datafile);
                numrpt = size(event_arr,2);
                valrptcnt = 0;
                maxwarp = 0;
                clear per clrchgind
                for rptlop = 1:numrpt
                    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) < 1010)) ~=0
                        if size(find(event_arr(:,rptlop) == 200)) ~=0
                            perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
                            perendind = find(event_arr(:,rptlop) == 24);
                            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                            begtimdum = time_arr(perbegind,rptlop)-100;
                            endtimdum = time_arr(perendind,rptlop);
                            if endtimdum > begtimdum
                                valrptcnt = valrptcnt + 1;
                                clrchgind(valrptcnt)=rptlop;
                                per(valrptcnt).begsmpind = begtimdum;
                                per(valrptcnt).endsmpind = endtimdum;
                                per(valrptcnt).begpos = 1;
                                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                                per(valrptcnt).allval = event_arr(:,rptlop);
                                per(valrptcnt).alltim = time_arr(:,rptlop);
                            end
                        end
                    end
                end
                clear cnd
                numrpt = size(per,2);
                for rptlop = 1:numrpt
                    cnd(rptlop)=per(rptlop).cnd;
                end
                evnnmb=2:2:size(eog_arr,1);
                oddnmb=1:2:size(eog_arr,1);
                clear x y
                cndlst=unique(cnd);
                for k=1:length(cndlst)
                    cndind=find(cnd==cndlst(k));
                    allind=clrchgind(cndind);
                    for l=1:length(allind)
                        x{k}(l)=mean(eog_arr(intersect(floor(((per(cndind(l)).begsmpind-1000)/samprate)*2):(floor((per(cndind(l)).endsmpind-1000)/samprate))*2,oddnmb),allind(l)));
                        y{k}(l)=mean(eog_arr(intersect(floor(((per(cndind(l)).begsmpind-1000)/samprate)*2):(floor((per(cndind(l)).endsmpind-1000)/samprate))*2,evnnmb),allind(l)));
                    end
                end
                %remove outlying points when calculating average eye position @ each location
                for k=1:length(x)
                    x{k}=x{k}(find(x{k}<mean(x{k}+std(x{k})) & x{k}>mean(x{k}-std(x{k}))));
                    y{k}=y{k}(find(x{k}<mean(x{k}+std(x{k})) & x{k}>mean(x{k}-std(x{k}))));
                    x{k}=x{k}(find(y{k}<mean(y{k}+std(y{k})) & y{k}>mean(y{k}-std(y{k}))));
                    y{k}=y{k}(find(y{k}<mean(y{k}+std(y{k})) & y{k}>mean(y{k}-std(y{k}))));
                end
                clear meanx meany
                for k=1:length(x)
                    meanx(k)=mean(x{k});
                end
                for k=1:length(y)
                    meany(k)=mean(y{k});
                end
                clear x y
                x=meanx; y=meany;
                meanxorigin = mean([x(6) x(2) x(1) x(5) x(9) ],2);
                xscale = mean([6/(x(8)-meanxorigin) 3/(x(4)-meanxorigin) 3/(abs(x(3)-meanxorigin)) 6/(abs(x(7)-meanxorigin))],2);
                meanyorigin = mean([y(7) y(3) y(1) y(4) y(8) ],2);
                yscale = mean([6/(y(6)-meanyorigin) 3/(y(2)-meanyorigin) 3/(abs(y(5)-meanyorigin)) 6/(abs(y(9)-meanyorigin))],2);
                numrpt = size(event_arr,2);
                valrptcnt = 0;
                clear per vpcind
                new_eog_arr=[];
                for rptlop = 1:numrpt
                    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) >= 1010)) ~=0
                        if size(find(event_arr(:,rptlop) == 200)) ~=0
                            perbegind = find(event_arr(:,rptlop) == 23,1,'first');
                            perendind = find(event_arr(:,rptlop) == 24,1,'first');
                            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                            begtimdum = time_arr(perbegind,rptlop);
                            endtimdum = time_arr(perendind,rptlop);
                            if endtimdum > begtimdum
                                valrptcnt = valrptcnt + 1;
                                vpcind(valrptcnt)=rptlop;
                                per(valrptcnt).begsmpind = begtimdum;
                                per(valrptcnt).endsmpind = endtimdum;
                                per(valrptcnt).begpos = 1;
                                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                                per(valrptcnt).allval = event_arr(:,rptlop);
                                per(valrptcnt).alltim = time_arr(:,rptlop);
                                new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                            end
                        end
                    end
                end
                eyedat = [];
                for trlop=1:size(per,2)
                    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
                    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
                    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
                    picstart=per(trlop).alltim(find(per(trlop).allval==23,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture start time relative to eye scan start
                    picend=per(trlop).alltim(find(per(trlop).allval==24,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture end time relative to eye scan start
                    
                    if picend > length(horeog)*5
                        picend =  length(horeog)*5;
                    end
                    eyedat{trlop}(1,:) = (horeog(ceil(picstart/5):floor(picend/5))) .* xscale;
                    eyedat{trlop}(2,:) = (vrteog(ceil(picstart/5):floor(picend/5))) .* yscale;
                end
                
                for i = 1:size(eyedat,2);
                    x = 24*eyedat{i}(1,:)+imageX/2;
                    y = 24*eyedat{i}(2,:)+imageY/2;
                    badx = find(x < -50 | x > imageX+50); %~1 dva leave margin of error
                    x(badx) = []; y(badx) = [];
                    bady = find(y < -50 | y > imageY+50); %~1 dva margin of error
                    x(bady) = []; y(bady) = [];
                    eyedat{i} = [x;y];
                end
            end
            
            fixationstats = ClusterFixation_Final(eyedat(1:2:end),samprate/1000);
            
            for ii = 1:length(fixationstats);
                if size(fixationstats{ii}.XY,2) >= 2000 %full viewing of image
                    img = zeros(imageY,imageX);
                    fixations = fixationstats{ii}.fixations;
                    if ~isempty(fixations)
                        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                            fixations(:,1) = [];
                        end
                        %remove fixations from outside image at least lose 1/2
                        %area these and upto almost 100% of area. Corrected in
                        %corrected batch for BCRW
                        badx = find((fixations(1,:) < 1)|(fixations(1,:) > imageX));
                        fixations(:,badx) = [];
                        bady = find((fixations(2,:) < 1)|(fixations(2,:) > imageY));
                        fixations(:,bady) = [];
                        for iii = 1:size(fixations,2)
                            x = fixations(1,iii);
                            y = fixations(2,iii);
                            C = sqrt((rr-x).^2+(cc-y).^2)<=IOR_area;
                            img(C) = 1;
                        end
                        if ~isempty(strfind(a.mat{aa},'MP'))
                            tag = 1;
                        elseif ~isempty(strfind(a.mat{aa},'TT'))
                            tag = 2;
                        elseif ~isempty(strfind(a.mat{aa},'JN'))
                            tag = 3;
                        elseif ~isempty(strfind(a.mat{aa},'IW'))
                            tag = 4;
                        end
                        percent_coverage(tag,count(tag)) = sum(sum(img))/numel(img);
                        count(tag) = count(tag)+1;
                    end
                end
            end
        end
    end
end
save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'Observed_Percent_Coverage_all'],'percent_coverage','count')
%% Calculate Percentage of Image Covered by BCRW
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;


IOR_dir = 'BCRW IOR TAU ';
IOR_taus = [0 50 35 25 17 12 7 3 1];

IOR_area = 48;
[rr,cc] = meshgrid(1:imageX,1:imageY);

for imset = 1:length(image_sets);
    for IOR = 1:length(IOR_taus);
        percent_coverage = NaN(100,144);
        count = 1;
        for i = 1:36;
            for t = 1:length(tags)
                load([image_dir '\' image_sets{imset} '\' IOR_dir num2str(IOR_taus(IOR)) '\' tags{t} '-' num2str(i) '-BCRW.mat'])
                disp(['Image set-' image_sets{imset} ' IOR_tau = ' num2str(IOR_taus(IOR)) ...
                    ' Image# ' num2str(i) ' Monkey ' tags{t}])
                for ii = 1:size(fixationtimes,1);
                    img = zeros(imageY,imageX);
                    tind = find(fixationtimes(ii,:,1) > 0);
                    for iii = 1:length(tind)
                        x = fixationtimes(ii,tind(iii),1);
                        y = fixationtimes(ii,tind(iii),2);
                        C = sqrt((rr-x).^2+(cc-y).^2)<=IOR_area;
                        img(C) = 1;
                    end
                    percent_coverage(ii,count) = sum(sum(img))/numel(img);
                end
                count = count+1;
            end
        end
        save([image_dir image_sets{imset} '\Percent_Image_Viewed-IOR_' num2str(IOR_taus(IOR)) '.mat'],'percent_coverage')
    end
end
%% Percent Coverage for Observed data from fixation files
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = ['MP';'TT';'JN';'IW'];
imageX = 800;
imageY = 600;
samprate = 5;

IOR_area = 48;
[rr,cc] = meshgrid(1:imageX,1:imageY);

count = ones(1,4);
percent_coverage = NaN(4,288);
for imset = 1:length(image_sets);
    cd([image_dir '\' image_sets{imset} '\'])
    a = what;
    for aa = 1:size(a.mat,1);
        if ~isempty(strfind(a.mat{aa},'fixation.mat'))
            disp(['Image set-' image_sets{imset} ' File ' a.mat{aa}])
            load(a.mat{aa})
            
            for ii = 1:2:length(fixationstats);
                img = zeros(imageY,imageX);
                fixations = fixationstats{ii}.fixations;
                if ~isempty(fixations)
                    if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                            fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                        fixations(:,1) = [];
                    end
                    %remove fixations from outside image at least lose 1/2
                    %area these and upto almost 100% of area. Corrected in
                    %corrected batch for BCRW
                    badx = find((fixations(1,:) < 1)|(fixations(1,:) > imageX));
                    fixations(:,badx) = [];
                    bady = find((fixations(2,:) < 1)|(fixations(2,:) > imageY));
                    fixations(:,bady) = [];
                    for iii = 1:size(fixations,2)
                        x = fixations(1,iii);
                        y = fixations(2,iii);
                        C = sqrt((rr-x).^2+(cc-y).^2)<=IOR_area;
                        img(C) = 1;
                    end
                    if ~isempty(strfind(a.mat{aa},'MP'))
                        tag = 1;
                    elseif ~isempty(strfind(a.mat{aa},'TT'))
                        tag = 2;
                    elseif ~isempty(strfind(a.mat{aa},'JN'))
                        tag = 3;
                    elseif ~isempty(strfind(a.mat{aa},'IW'))
                        tag = 4;
                    end
                    percent_coverage(tag,count(tag)) = sum(sum(img))/numel(img);
                    count(tag) = count(tag)+1;
                end
            end
        end
    end
end
save(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'Observed_Percent_Coverage_10s'],'percent_coverage','count')


%% Calculate Percentage of Image Covered by BCRW while limiting number of
% fixations to bserved number for each image
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;


IOR_dir = 'BCRW IOR TAU ';
IOR_taus = [0 50 35 25 17 12 7 3 1];

IOR_area = 48;
[rr,cc] = meshgrid(1:imageX,1:imageY);

for imset = 1:length(image_sets);
    cd([image_dir image_sets{imset}])
    a = what;
    fixationfiles = {};
    for aa = 1:length(a.mat);
        if ~isempty(strfind(a.mat{aa},'-fixation.mat'))
            if ~isempty(strfind(a.mat{aa},'MP'))
                fixationfiles{1} = a.mat{aa};
            elseif ~isempty(strfind(a.mat{aa},'TT'))
                fixationfiles{2} = a.mat{aa};
            elseif ~isempty(strfind(a.mat{aa},'JN'))
                fixationfiles{3} = a.mat{aa};
            elseif ~isempty(strfind(a.mat{aa},'IW'))
                fixationfiles{4} = a.mat{aa};
            end
        end
    end
    for IOR = 1:length(IOR_taus);
        percent_coverage = NaN(100,144);
        count = 1;
        for i = 1:36;
            for t = 1:length(tags)
                load([image_dir '\' image_sets{imset} '\' IOR_dir num2str(IOR_taus(IOR)) '\' tags{t} '-' num2str(i) '-BCRW.mat'])
                disp(['Image set-' image_sets{imset} ' IOR_tau = ' num2str(IOR_taus(IOR)) ...
                    ' Image# ' num2str(i) ' Monkey ' tags{t}])
                load(fixationfiles{t})
                for ii = 1:size(fixationtimes,1);
                    img = zeros(imageY,imageX);
                    tind = find(fixationtimes(ii,:,1) > 0);
                    
                    fixations = fixationstats{2*i-1}.fixations;
                    if isempty(fixations)
                        continue
                    end
                    if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                            fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                        fixations(:,1) = [];
                    end
                    if isempty(fixations)
                        continue
                    end
                    %remove fixations from outside image at least lose 1/2
                    %area these and upto almost 100% of area. Corrected in
                    %corrected batch for BCRW
                    badx = find((fixations(1,:) < 1)|(fixations(1,:) > imageX));
                    fixations(:,badx) = [];
                    bady = find((fixations(2,:) < 1)|(fixations(2,:) > imageY));
                    fixations(:,bady) = [];
                    if length(tind) > size(fixations,2)
                        tind = tind(1:size(fixations,2));
                    end
                    for iii = 1:length(tind)
                        x = fixationtimes(ii,tind(iii),1);
                        y = fixationtimes(ii,tind(iii),2);
                        C = sqrt((rr-x).^2+(cc-y).^2)<=IOR_area;
                        img(C) = 1;
                    end
                    percent_coverage(ii,count) = sum(sum(img))/numel(img);
                end
                count = count+1;
            end
        end
        save([image_dir image_sets{imset} '\Percent_Image_Viewed_corrected-IOR_' num2str(IOR_taus(IOR)) '.mat'],'percent_coverage')
    end
end
%% Plot and analyze percent coverage for various IOR values in BCRW
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};

IOR_taus = [0 50 35 25 17 12 7 3 1];
IORTau =  [0.01 1./[50 35 25 17 12 7 3 1]];%corrected for plotting

labels ={};
for ir = 1:length(IOR_taus);
    if ir == 1
        labels{ir} = '0';
    elseif ir == 9;
        labels{ir} = '1';
    else
        labels{ir} = ['1/' num2str(1/IORTau(ir))];
    end
end

all_pc = cell(1,length(IOR_taus));
for imset = 1:length(image_sets)
    for ir = 1:length(IOR_taus);
        filename = [image_dir image_sets{imset} '\Percent_Image_Viewed_corrected-IOR_' num2str(IOR_taus(ir)) '.mat'];
        if exist(filename)
            load(filename);
            pc = percent_coverage(:,1:144); %some are sized incorrectly but full of NaNs.
            all_pc{ir} = [all_pc{ir}; nanmean(pc)];
        end
    end
end

mean_pcs = [];
std_pcs = [];
num_pcs = [];
for ir = 1:length(IOR_taus);
    temp = [];
    for t = 1:4;
        temp =  [temp; all_pc{ir}(t:4:end)];
    end
    temp = mean(temp);
    num_pcs(ir) = sqrt(sum(sum(~isnan(temp)))/4);
    mean_pcs(ir) = nanmean(temp(1:end));
    std_pcs(ir) = nanstd(temp(1:end));
end
figure
%  hold on
errorbar(log(IORTau),100*mean_pcs,100*std_pcs./sqrt(num_pcs),'b')
xlabel('Log(IOR_tau)')
ylabel('Percent Coverage')
set(gca,'Xtick',[log(IORTau)])
set(gca,'XtickLabel',labels)
%% Compare observed to BCRW data
load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\Observed_Percent_Coverage_10s.mat'])

p_vals = [];
percent_coverage = mean(percent_coverage);
for ir = 1:length(all_pc);
    ap = [];
    for t = 1:length(tags);
        ap = [ap; all_pc{ir}(t:4:end)];
    end
    ap = mean(ap);
    [~,p] = ttest2(percent_coverage,ap);
    p_vals(ir) = p;
end
%% Correlated observed to BCRW data
tags = {'MP','TT','JN','IW'};

image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;


IOR_taus = [0 50 35 25 17 12 7 3 1];
bcrw_data = cell(length(IOR_taus),4);

for imset = 1:length(image_sets);
    cd([image_dir image_sets{imset}])
    for ir = 1:length(IOR_taus)
        load(['Percent_Image_Viewed_corrected-IOR_' num2str(IOR_taus(ir)) '.mat'])
        pc = nanmean(percent_coverage);
        for i = 1:4;
            bcrw_data{ir,i} = [bcrw_data{ir,i} pc(i:4:end)];
        end
    end
end

load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\Observed_Percent_Coverage_10s.mat'])
percent_coverage = percent_coverage(:,1:2:end);

%% CURVE FITTING OMG NOOOOOOOO BUT JUST GOING TO SEE INCASE
load(['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\Observed_Percent_Coverage_10s.mat'])

x = log(1./IOR_taus(4:end));
y = 100*mean_pcs(4:end);
p = polyfit(x,y,1);

c = corrcoef(x,y);
c=c(2);
estimate_IOR= 1/exp((nanmean(nanmean(100*percent_coverage))-p(2))/p(1))
%% Histograms of observed vs Predicted Coverage
for ir = 1:length(IOR_taus);
    figure
    for i = 1:4
        subplot(2,2,i)
        p = percent_coverage(i,:);
        b = bcrw_data{ir,i};
        r = corrcoef(p,b);
        r = r(2);
        pf = polyfit(p,b,1);
        plot(p,b,'.');
        title([tags{i} ' corrcoef: ' num2str(r) 'slope: ' num2str(pf(1))])
        xlabel('Observed')
        ylabel('Predicted')
        axis square
    end
    subtitle(['IOR_{tau} ' num2str(IOR_taus(ir))])
end

for ir = 1:length(IOR_taus);
    figure
    for i = 1:4
        p = percent_coverage(i,:);
        b = bcrw_data{ir,i};
        subplot(2,4,i)
        hist(p,20)
        xlim([0 0.5])
        xlabel('Observed')
        title(['Mean=' num2str(nanmean(p)) ' std=' num2str(nanstd(p))]);
        subplot(2,4,4+i)
        hist(b,20)
        xlim([0 0.5])
        xlabel('Predicted')
        title(['Mean=' num2str(nanmean(b)) ' std=' num2str(nanstd(b))]);
    end
    subtitle(['IOR_{tau} ' num2str(IOR_taus(ir))])
end