function SalienceIOR(FIXATIONFILE,distancethreshold,imageX,imageY)
% created by Seth Koenig 11/21/2012

% function determines the rate of return fixations, the time between return
% fixations, time within trial of return, and salience at returned
% location.

% Inputs:
%   FIXATIONFILE: Fixations extracted from cortex e.g. MP120606_1-fixation.mat
%   ImageX,ImageY: x and y dimensions of images
%   Pairings: take all pairing or only closest unique pairing between first
%   and second fixation


%Outputs:
%   A .mat file named [FIXATIONFILE(1:end-13) '-SalienceIOR']
%   containg saccade statistics. See variable statvariablenames for
%   detailed explanation of variables in .mat file.

if nargin < 1
    error(['Not enough inputs: function requires FixationFile,'...
        'distance threhsold,imageX, imageY, and pairings.'])
end
if nargin < 2
    distancethreshold = [0  24 48 72 96  120  144 168 200 400 800;...%in pixels 24 pixels/dva
        24 48 72 96 120 144  168 200 400 800 10000];%in pixels 24 pixels/dva
end
if nargin < 4
    imageX = 800;
    imageY = 600;
end

load(FIXATIONFILE);
matfiles = what;
saliencemapfiles = [NaN;NaN];
for i = 1:length(matfiles.mat);
    str = strfind(matfiles.mat{i},'saliencemap.mat');
    if ~isempty(str)
        dash = strfind(matfiles.mat{i},'-');
        if ~isempty(str2num(matfiles.mat{i}(dash(1)-1)))
            saliencemapfiles = [saliencemapfiles [i;str2num(matfiles.mat{i}(1:dash(1)-1))]];
        end
    end
end
saliencemapfiles(:,1) = [];
[~,si] = sort(saliencemapfiles(2,:));
saliencemapfiles = saliencemapfiles(1,si);

returnfixsal = cell(1,size(distancethreshold,2));
%fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist fix1dur fix2dur fixnum1 fixnum2
count = ones(1,size(distancethreshold,2));
for i = 1:size(distancethreshold,2)
    returnfixsal{i} = NaN(250,13);
end

for cndlop=1:2:length(fixationstats)
    reindexed = (cndlop+1)/2;
    load(matfiles.mat{saliencemapfiles(reindexed)},'fullmap');
    saliencemap = fullmap;
    fixations = fixationstats{cndlop}.fixations;
    if ~isempty(fixations)
        fixationtimes = fixationstats{cndlop}.fixationtimes;
        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
            fixations(:,1) = [];
            fixationtimes(:,1) = [];
        end
        
        N=size(fixations,2);
        [x,y]=meshgrid(1:N);
        i=find(ones(N)-eye(N)); %forms pairs except for self-pairing
        i=[x(i), y(i)];
        i(i(:,1) > i(:,2),:) = []; %repeat pairs
        i(i(:,1)+1 == i(:,2),:) = []; %removes consecutive in time pairs
        dist =sqrt((fixations(1,i(:,1))-fixations(1,i(:,2))).^2 +...
            (fixations(2,i(:,1))-fixations(2,i(:,2))).^2);
        wcount = 1;
        pairs = NaN(ceil(size(fixations,2)/2),3);
        while ~isempty(i);
            [minn,mind] = min(dist);
            minn = minn(1); %if multiple just take the first
            mind = mind(1); %if multiple just take the first
            
            middlefixes = i(mind,1)+1:i(mind,2)-1;
            middist = sqrt((fixations(1,i(mind,1))-fixations(1,middlefixes)).^2+...
                (fixations(2,i(mind,2))-fixations(2,middlefixes)).^2);
            
            if any(middist >= 10*24) %so at least 1 10 dva saccade away from area
                tind = find(minn > distancethreshold(1,:) & minn <= distancethreshold(2,:));
                if ~isempty(tind);
                    pairs(wcount,:) = [i(mind,:) tind];
                end
                %remove all instances of fixations with prior or return number
                [rmvind1,~] = find(i(:,1) == i(mind,1));
                [rmvind2,~] = find(i(:,2) == i(mind,1));
                [rmvind3,~] = find(i(:,1) == i(mind,2));
                [rmvind4,~] = find(i(:,2) == i(mind,2));
                rmvind = [rmvind1; rmvind2; rmvind3; rmvind4];
                rmvind = unique(rmvind);
                i(rmvind,:) = [];
                dist(rmvind) = [];
                wcount = wcount+1;
            else %then remove pair
                dist(mind) = [];
                i(mind,:) = [];
            end
        end
        pairs(isnan(pairs(:,1)),:) = [];
        
        for i = 1:size(pairs,1);
            
            spot = [ceil(fixations(:,pairs(i,1))) ceil(fixations(:,pairs(i,2)))];
            spot(2,:) = imageY-spot(2,:); %location
            spott = [(fixationtimes(1,pairs(i,1))+fixationtimes(2,pairs(i,1)))/2 ...
                (fixationtimes(1,pairs(i,2))+fixationtimes(2,pairs(i,2)))/2]; %time 
            dist = sqrt((spot(1,1)-spot(1,2))^2+(spot(2,1)-spot(2,2))^2);
            spot(spot < 1) = 1;
            spot(1,spot(1,:) > imageX) = imageX;
            spot(2,spot(2,:) > imageY) = imageY;
            
            returnfixsal{pairs(i,3)}(count(pairs(i,3)),:) = [...
                spot(1,1) spot(2,1) spott(1) saliencemap(spot(2,1),spot(1,1))...
                spot(1,2) spot(2,2) spott(2) saliencemap(spot(2,2),spot(1,2))...
                dist diff(fixationtimes(:,pairs(i,1)))+1 ...
                diff(fixationtimes(:,pairs(i,2)))+1 pairs(i,1) pairs(i,2)];
            %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist fix1dur fix2dur fixnum1 fixnum2
            
            count(pairs(i,3)) = count(pairs(i,3))+1;
        end
    end
end
for i = 1:size(distancethreshold,2)
    returnfixsal{i}(isnan(returnfixsal{i}(:,1)),:) = [];
end

IORvariablenames = {
    'returnfixsal: [  %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal...';
    'fixdist fix1dur fix2dur fixnum1 fixnum2]';
    };

save([FIXATIONFILE(1:end-13) '-SalienceIOR'],'returnfixsal','IORvariablenames')
end