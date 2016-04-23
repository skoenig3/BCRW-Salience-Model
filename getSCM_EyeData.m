function [eyedat,trialtype]  = getSCM_EyeData(datafile,samprate)
%updated from getEyeData.m by Seth Konig 5/20/15
%function imports SCM eye data
%updates mainly improvement in calibration function using cp2tform instead
%of determining the scale linearly. Tries to get best calibration function.
%Other update is to correct for block error during Irwin's set 6 (skipped
%35th 2nd presentation on break fixation errro), and attemtp to catch any
%similar errors.

samprate = samprate*1000;
%sub function convert raw eye tracking data into x & y coordinates
[time_arr,event_arr,eog_arr,~,~,~]  = get_ALLdata(datafile);
numrpt = size(event_arr,2);
valrptcnt = 0;
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
%clear cnd
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

% old calibration code--removed 5/20/15
% clear x y
% x=meanx; y=meany;
% meanxorigin = mean([x(6) x(2) x(1) x(5) x(9) ],2);
% xscale = mean([6/(x(8)-meanxorigin) 3/(x(4)-meanxorigin) 3/(abs(x(3)-meanxorigin)) 6/(abs(x(7)-meanxorigin))],2);
% meanyorigin = mean([y(7) y(3) y(1) y(4) y(8) ],2);
% yscale = mean([6/(y(6)-meanyorigin) 3/(y(2)-meanyorigin) 3/(abs(y(5)-meanyorigin)) 6/(abs(y(9)-meanyorigin))],2);

%new calibration code
controlx = [0 0 -3 3 0 0 -6 6 0];
controly = [0 3 0 0 -3 6 0 0 -6];
tform = get_calibration_fcn([controlx; controly],[meanx;meany]);


numrpt = size(event_arr,2);
valrptcnt = 0;
clear per vpcind
new_eog_arr=[];
trialtype = '';
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
                cnd = find(event_arr(:,rptlop) >= 1010 & event_arr(:,rptlop) < 2000);
                type = event_arr(cnd-1,rptlop);
                if type == 1;
                    trialtype = [trialtype 'n']; %novel
                elseif type == 2
                    trialtype = [trialtype 'f']; %familiar/repeat
                elseif type == 3
                    trialtype = [trialtype 'r']; %object replaced
                elseif type == 4;
                    trialtype = [trialtype 'm']; %object moved
                else
                    trialtype = [trialtype 'e']; %some kind of error occured
                    disp('e')
                    disp(type)
                end
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
    eyedat{trlop}(1,:) = (horeog(ceil(picstart/5):floor(picend/5)));% .* xscale; %removed 5/20/15
    eyedat{trlop}(2,:) = (vrteog(ceil(picstart/5):floor(picend/5)));% .* yscale; %removed 5/20/15
end

%Block error with Irwin's file skipped 35th familiar image, was supposed to
%be moved
if length(eyedat) < 72
    if strcmpi(datafile,'IW100707.1')
        trialtype = [trialtype(1:69) 'e' trialtype(70:71)];
        eyedat = [eyedat(1:69) {NaN(2,1)} eyedat(70:71)];
    else
        error(['Not all images were displayed or other error occurred. File: ' datafile])
    end
end

%---Recalibrate and automatically scale eye data---%
% added 5/20/15
for eye = 1:length(eyedat)
    if all(~isnan(eyedat{eye}))
        x = eyedat{eye}(1,:);
        y = eyedat{eye}(2,:);
        [x,y] = tformfwd(tform,x,y);
        eyedat{eye} = [x;y];
    end
end

%remove eye data outside of image
for i = 1:size(eyedat,2);
    x = eyedat{i}(1,:)*24+400;
    y = eyedat{i}(2,:)*24+300;
    %code shouldn't be necessary since image start is already time 0
    %     tstart = find(x > 300 & x < 500 & y > 200 & y < 400);
    %     tstart = tstart(1);
    %     x = x(tstart:end); y = y(tstart:end);
    if length(x) > 2000; %1st 10 seconds
        x(2001:end) = []; y(2001:end) = [];
    end
    badx = find(x < -25 | x > 825); %~1 dva leave margin of error, was 2 dva SDK 5/26 but calibration looks better now
    x(badx) = []; y(badx) = [];
    bady = find(y < -25 | y > 625); %~1 dva margin of error, was 2 dva SDK 5/26 but calibration looks better now
    x(bady) = []; y(bady) = [];
    x = (x-400)/24; y = (y-300)/24;
    eyedat{i} = [x;y];
end
end