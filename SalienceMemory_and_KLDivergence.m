% Test Salience and KL divergence between monkeys that view novel vs
% repeated images. Currently still incomplete.
%---[9] Combine BCRW and Calculate Goodness of Fit-KL Divergence---%
scm_image_dir = 'C:\Users\skoenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','JN','IW'};
imageX = 800;
imageY = 600;
binsize = 25; %24 pixels per dva but using 25 cuz divides nicely into 600 & 800
f = fspecial('gaussian',[256,256],24);
KLnorm = NaN(36*length(image_sets),6,length(tags));
KLshuff = NaN(36*length(image_sets),6,length(tags));
novsalfix = NaN(36*length(image_sets),75,length(tags));
repsalfix = NaN(36*length(image_sets),75,length(tags));

prebiasPDF15 = zeros(imageY,imageX);
prebiasPDF610 = zeros(imageY,imageX);
prebiasPDF1115 = zeros(imageY,imageX);
prebiasPDF1620 = zeros(imageY,imageX);
prebiasPDF2125 = zeros(imageY,imageX);
prebiasPDF = zeros(imageY,imageX);

postbiasPDF15 = zeros(imageY,imageX);
postbiasPDF610 = zeros(imageY,imageX);
postbiasPDF1115 = zeros(imageY,imageX);
postbiasPDF1620 = zeros(imageY,imageX);
postbiasPDF2125 = zeros(imageY,imageX);
postbiasPDF = zeros(imageY,imageX);

prebiasshuffPDF15 = zeros(imageY,imageX);
prebiasshuffPDF610 = zeros(imageY,imageX);
prebiasshuffPDF1115 = zeros(imageY,imageX);
prebiasshuffPDF1620 = zeros(imageY,imageX);
prebiasshuffPDF2125 = zeros(imageY,imageX);
prebiasshuffPDF = zeros(imageY,imageX);

postbiasshuffPDF15 = zeros(imageY,imageX);
postbiasshuffPDF610 = zeros(imageY,imageX);
postbiasshuffPDF1115 = zeros(imageY,imageX);
postbiasshuffPDF1620 = zeros(imageY,imageX);
postbiasshuffPDF2125 = zeros(imageY,imageX);
postbiasshuffPDF = zeros(imageY,imageX);

for imset = 1:length(image_sets);
    disp(['Image set-' num2str(image_sets{imset})])
    dirName = [scm_image_dir image_sets{imset}];
    cd(dirName)
    matfiles = what;
    
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
    
    for i = 1:36;
        allfixations = zeros(imageY,imageX);
        allBCRW = zeros(imageY,imageX);
        index = (imset-1)*36+i;
        
        load(matfiles.mat{saliencemapfiles(i)})
        
        for t = 1:length(tags)
            if eyedatafiles(t) ~= 0;
                
                load(matfiles.mat{eyedatafiles(t)})
                if length(trialtype) >= 2*i
                    if strcmp(trialtype(2*i),'f')
                        fixations = fixationstats{i*2}.fixations;
                        fixationtimes = fixationstats{i*2}.fixationtimes;
                        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                            fixations(:,1) = [];
                            fixationtimes(:,1) = [];
                        end
                        if length(fixations) >= 5;
                            fixations = fixationstats{i*2-1}.fixations;
                            fixationtimes = fixationstats{i*2-1}.fixationtimes;
                            if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                                    fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                                fixations(:,1) = [];
                                fixationtimes(:,1) = [];
                            end
                            prefixPDF15 = zeros(imageY,imageX);
                            prefixPDF610 = zeros(imageY,imageX);
                            prefixPDF1115 = zeros(imageY,imageX);
                            prefixPDF1620 = zeros(imageY,imageX);
                            prefixPDF2125 = zeros(imageY,imageX);
                            
                            prefixshuffPDF15 = zeros(imageY,imageX);
                            prefixshuffPDF610 = zeros(imageY,imageX);
                            prefixshuffPDF1115 = zeros(imageY,imageX);
                            prefixshuffPDF1620 = zeros(imageY,imageX);
                            prefixshuffPDF2125 = zeros(imageY,imageX);
                            for ii = 1:length(fixations)
                                spot = ceil(fixations(:,ii));
                                spot(2) = imageY-spot(2);
                                spot(spot < 1) = 1;
                                spot(1,spot(1) > imageX) = imageX;
                                spot(2,spot(2) > imageY) = imageY;
                                novsalfix(index,ii,t) = fullmap(spot(2),spot(1));
                                randspot = [ceil(rand*800) ceil(rand*600)];
                                if ii <= 5
                                    prefixPDF15(spot(2),spot(1)) = prefixPDF15(spot(2),spot(1)) + 1;
                                elseif ii <= 10
                                    prefixPDF610(spot(2),spot(1)) = prefixPDF610(spot(2),spot(1)) + 1;
                                elseif ii <= 15
                                    prefixPDF1115(spot(2),spot(1)) = prefixPDF1115(spot(2),spot(1)) + 1;
                                elseif ii <= 20
                                    prefixPDF1620(spot(2),spot(1)) = prefixPDF1620(spot(2),spot(1)) + 1;
                                elseif ii <= 25
                                    prefixPDF2125(spot(2),spot(1)) = prefixPDF2125(spot(2),spot(1)) + 1;
                                end
                                if ii <= 5
                                    prefixshuffPDF15(randspot(2),randspot(1)) = prefixshuffPDF15(randspot(2),randspot(1)) + 1;
                                elseif ii <= 10
                                    prefixshuffPDF610(randspot(2),randspot(1)) = prefixshuffPDF610(randspot(2),randspot(1)) + 1;
                                elseif ii <= 15
                                    prefixshuffPDF1115(randspot(2),randspot(1)) = prefixshuffPDF1115(randspot(2),randspot(1)) + 1;
                                elseif ii <= 20
                                    prefixshuffPDF1620(randspot(2),randspot(1)) = prefixshuffPDF1620(randspot(2),randspot(1)) + 1;
                                elseif ii <= 25
                                    prefixshuffPDF2125(randspot(2),randspot(1)) = prefixshuffPDF2125(randspot(2),randspot(1)) + 1;
                                end
                            end
                            prefixPDF = prefixPDF15 + prefixPDF610 + prefixPDF1115 + ...
                                prefixPDF1620 + prefixPDF2125;
                            prefixshuffPDF = prefixshuffPDF15 + prefixshuffPDF610 + prefixshuffPDF1115 + ...
                                prefixshuffPDF1620 + prefixshuffPDF2125;
                            fixations = fixationstats{i*2}.fixations;
                            fixationtimes = fixationstats{i*2}.fixationtimes;
                            if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                                    fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                                fixations(:,1) = [];
                                fixationtimes(:,1) = [];
                            end
                            
                            postfixPDF15 = zeros(imageY,imageX);
                            postfixPDF610 = zeros(imageY,imageX);
                            postfixPDF1115 = zeros(imageY,imageX);
                            postfixPDF1620 = zeros(imageY,imageX);
                            postfixPDF2125 = zeros(imageY,imageX);
                            
                            postfixshuffPDF15 = zeros(imageY,imageX);
                            postfixshuffPDF610 = zeros(imageY,imageX);
                            postfixshuffPDF1115 = zeros(imageY,imageX);
                            postfixshuffPDF1620 = zeros(imageY,imageX);
                            postfixshuffPDF2125 = zeros(imageY,imageX);
                            for ii = 1:length(fixations)
                                spot = ceil(fixations(:,ii));
                                spot(2) = imageY-spot(2);
                                spot(spot < 1) = 1;
                                spot(1,spot(1) > imageX) = imageX;
                                spot(2,spot(2) > imageY) = imageY;
                                repsalfix(index,ii,t) = fullmap(spot(2),spot(1));
                                randspot = [ceil(rand*800) ceil(rand*600)];
                                if ii <= 5
                                    postfixPDF15(spot(2),spot(1)) = postfixPDF15(spot(2),spot(1)) + 1;
                                elseif ii <= 10
                                    postfixPDF610(spot(2),spot(1)) = postfixPDF610(spot(2),spot(1)) + 1;
                                elseif ii <= 15
                                    postfixPDF1115(spot(2),spot(1)) = postfixPDF1115(spot(2),spot(1)) + 1;
                                elseif ii <= 20
                                    postfixPDF1620(spot(2),spot(1)) = postfixPDF1620(spot(2),spot(1)) + 1;
                                elseif ii <= 25
                                    postfixPDF2125(spot(2),spot(1)) = postfixPDF2125(spot(2),spot(1)) + 1;
                                end
                                if ii <= 5
                                    postfixshuffPDF15(randspot(2),randspot(1)) = postfixshuffPDF15(randspot(2),randspot(1)) + 1;
                                elseif ii <= 10
                                    postfixshuffPDF610(randspot(2),randspot(1)) = postfixshuffPDF610(randspot(2),randspot(1)) + 1;
                                elseif ii <= 15
                                    postfixshuffPDF1115(randspot(2),randspot(1)) = postfixshuffPDF1115(randspot(2),randspot(1)) + 1;
                                elseif ii <= 20
                                    postfixshuffPDF1620(randspot(2),randspot(1)) = postfixshuffPDF1620(randspot(2),randspot(1)) + 1;
                                elseif ii <= 25
                                    postfixshuffPDF2125(randspot(2),randspot(1)) = postfixshuffPDF2125(randspot(2),randspot(1)) + 1;
                                end
                            end
                            postfixPDF = postfixPDF15 + postfixPDF610 + postfixPDF1115 + ...
                                postfixPDF1620 + postfixPDF2125;
                            postfixshuffPDF = postfixshuffPDF15 + postfixshuffPDF610 + postfixshuffPDF1115 + ...
                                postfixshuffPDF1620 + postfixshuffPDF2125;
                            
                            prebiasPDF15 = prebiasPDF15+prefixPDF15;
                            prebiasPDF610 = prebiasPDF610+prefixPDF610;
                            prebiasPDF1115 = prebiasPDF1115+prefixPDF1115;
                            prebiasPDF1620 = prebiasPDF1620+prefixPDF1620;
                            prebiasPDF2125 = prebiasPDF2125+prefixPDF2125;
                            prebiasPDF = prebiasPDF+prefixPDF;
                            
                            postbiasPDF15 = postbiasPDF15+postfixPDF15;
                            postbiasPDF610 = postbiasPDF610+postfixPDF610;
                            postbiasPDF1115 = postbiasPDF1115+postfixPDF1115;
                            postbiasPDF1620 = postbiasPDF1620+postfixPDF1620;
                            postbiasPDF2125 = postbiasPDF2125+postfixPDF2125;
                            postbiasPDF = postbiasPDF+postfixPDF;
                            
                            prefixPDF15 = imfilter(prefixPDF15,f);
                            prefixPDF610 = imfilter(prefixPDF610,f);
                            prefixPDF1115 = imfilter(prefixPDF1115,f);
                            prefixPDF1620 = imfilter(prefixPDF1620,f);
                            prefixPDF2125 = imfilter(prefixPDF2125,f);
                            prefixPDF = imfilter(prefixPDF,f);
                            
                            prefixPDF15 = bin2(prefixPDF15,binsize,binsize);
                            prefixPDF610 = bin2(prefixPDF610,binsize,binsize);
                            prefixPDF1115 = bin2(prefixPDF1115,binsize,binsize);
                            prefixPDF1620 = bin2(prefixPDF1620,binsize,binsize);
                            prefixPDF2125 = bin2(prefixPDF2125,binsize,binsize);
                            prefixPDF = bin2(prefixPDF,binsize,binsize);
                            
                            prefixPDF15(prefixPDF15 == 0) = eps;
                            prefixPDF610(prefixPDF610 == 0) = eps;
                            prefixPDF1115(prefixPDF1115 == 0) = eps;
                            prefixPDF1620(prefixPDF1620 == 0) = eps;
                            prefixPDF2125(prefixPDF2125 == 0) = eps;
                            prefixPDF(prefixPDF == 0) = eps;
                            
                            prefixPDF15 = prefixPDF15/sum(sum(prefixPDF15));
                            prefixPDF610 = prefixPDF610/sum(sum(prefixPDF610));
                            prefixPDF1115 = prefixPDF1115/sum(sum(prefixPDF1115));
                            prefixPDF1620 = prefixPDF1620/sum(sum(prefixPDF1620));
                            prefixPDF2125 = prefixPDF2125/sum(sum(prefixPDF2125));
                            prefixPDF = prefixPDF/sum(sum(prefixPDF));
                            
                            postfixPDF15 = imfilter(postfixPDF15,f);
                            postfixPDF610 = imfilter(postfixPDF610,f);
                            postfixPDF1115 = imfilter(postfixPDF1115,f);
                            postfixPDF1620 = imfilter(postfixPDF1620,f);
                            postfixPDF2125 = imfilter(postfixPDF2125,f);
                            postfixPDF = imfilter(postfixPDF,f);
                            
                            postfixPDF15 = bin2(postfixPDF15,binsize,binsize);
                            postfixPDF610 = bin2(postfixPDF610,binsize,binsize);
                            postfixPDF1115 = bin2(postfixPDF1115,binsize,binsize);
                            postfixPDF1620 = bin2(postfixPDF1620,binsize,binsize);
                            postfixPDF2125 = bin2(postfixPDF2125,binsize,binsize);
                            postfixPDF = bin2(postfixPDF,binsize,binsize);
                            
                            postfixPDF15(postfixPDF15 == 0) = eps;
                            postfixPDF610(postfixPDF610 == 0) = eps;
                            postfixPDF1115(postfixPDF1115 == 0) = eps;
                            postfixPDF1620(postfixPDF1620 == 0) = eps;
                            postfixPDF2125(postfixPDF2125 == 0) = eps;
                            postfixPDF(postfixPDF == 0) = eps;
                            
                            postfixPDF15 = postfixPDF15/sum(sum(postfixPDF15));
                            postfixPDF610 = postfixPDF610/sum(sum(postfixPDF610));
                            postfixPDF1115 = postfixPDF1115/sum(sum(postfixPDF1115));
                            postfixPDF1620 = postfixPDF1620/sum(sum(postfixPDF1620));
                            postfixPDF2125 = postfixPDF2125/sum(sum(postfixPDF2125));
                            postfixPDF = postfixPDF/sum(sum(postfixPDF));
                            
                            
                            prebiasshuffPDF15 = prebiasshuffPDF15+prefixshuffPDF15;
                            prebiasshuffPDF610 = prebiasshuffPDF610+prefixshuffPDF610;
                            prebiasshuffPDF1115 = prebiasshuffPDF1115+prefixshuffPDF1115;
                            prebiasshuffPDF1620 = prebiasshuffPDF1620+prefixshuffPDF1620;
                            prebiasshuffPDF2125 = prebiasshuffPDF2125+prefixshuffPDF2125;
                            prebiasshuffPDF = prebiasshuffPDF+prefixshuffPDF;
                            
                            postbiasshuffPDF15 = postbiasshuffPDF15+postfixshuffPDF15;
                            postbiasshuffPDF610 = postbiasshuffPDF610+postfixshuffPDF610;
                            postbiasshuffPDF1115 = postbiasshuffPDF1115+postfixshuffPDF1115;
                            postbiasshuffPDF1620 = postbiasshuffPDF1620+postfixshuffPDF1620;
                            postbiasshuffPDF2125 = postbiasshuffPDF2125+postfixshuffPDF2125;
                            postbiasshuffPDF = postbiasshuffPDF+postfixshuffPDF;
                            
                            prefixshuffPDF15 = imfilter(prefixshuffPDF15,f);
                            prefixshuffPDF610 = imfilter(prefixshuffPDF610,f);
                            prefixshuffPDF1115 = imfilter(prefixshuffPDF1115,f);
                            prefixshuffPDF1620 = imfilter(prefixshuffPDF1620,f);
                            prefixshuffPDF2125 = imfilter(prefixshuffPDF2125,f);
                            prefixshuffPDF = imfilter(prefixshuffPDF,f);
                            
                            prefixshuffPDF15 = bin2(prefixshuffPDF15,binsize,binsize);
                            prefixshuffPDF610 = bin2(prefixshuffPDF610,binsize,binsize);
                            prefixshuffPDF1115 = bin2(prefixshuffPDF1115,binsize,binsize);
                            prefixshuffPDF1620 = bin2(prefixshuffPDF1620,binsize,binsize);
                            prefixshuffPDF2125 = bin2(prefixshuffPDF2125,binsize,binsize);
                            prefixshuffPDF = bin2(prefixshuffPDF,binsize,binsize);
                            
                            prefixshuffPDF15(prefixshuffPDF15 == 0) = eps;
                            prefixshuffPDF610(prefixshuffPDF610 == 0) = eps;
                            prefixshuffPDF1115(prefixshuffPDF1115 == 0) = eps;
                            prefixshuffPDF1620(prefixshuffPDF1620 == 0) = eps;
                            prefixshuffPDF2125(prefixshuffPDF2125 == 0) = eps;
                            prefixshuffPDF(prefixshuffPDF == 0) = eps;
                            
                            prefixshuffPDF15 = prefixshuffPDF15/sum(sum(prefixshuffPDF15));
                            prefixshuffPDF610 = prefixshuffPDF610/sum(sum(prefixshuffPDF610));
                            prefixshuffPDF1115 = prefixshuffPDF1115/sum(sum(prefixshuffPDF1115));
                            prefixshuffPDF1620 = prefixshuffPDF1620/sum(sum(prefixshuffPDF1620));
                            prefixshuffPDF2125 = prefixshuffPDF2125/sum(sum(prefixshuffPDF2125));
                            prefixshuffPDF = prefixshuffPDF/sum(sum(prefixshuffPDF));
                            
                            postfixshuffPDF15 = imfilter(postfixshuffPDF15,f);
                            postfixshuffPDF610 = imfilter(postfixshuffPDF610,f);
                            postfixshuffPDF1115 = imfilter(postfixshuffPDF1115,f);
                            postfixshuffPDF1620 = imfilter(postfixshuffPDF1620,f);
                            postfixshuffPDF2125 = imfilter(postfixshuffPDF2125,f);
                            postfixshuffPDF = imfilter(postfixshuffPDF,f);
                            
                            postfixshuffPDF15 = bin2(postfixshuffPDF15,binsize,binsize);
                            postfixshuffPDF610 = bin2(postfixshuffPDF610,binsize,binsize);
                            postfixshuffPDF1115 = bin2(postfixshuffPDF1115,binsize,binsize);
                            postfixshuffPDF1620 = bin2(postfixshuffPDF1620,binsize,binsize);
                            postfixshuffPDF2125 = bin2(postfixshuffPDF2125,binsize,binsize);
                            postfixshuffPDF = bin2(postfixshuffPDF,binsize,binsize);
                            
                            postfixshuffPDF15(postfixshuffPDF15 == 0) = eps;
                            postfixshuffPDF610(postfixshuffPDF610 == 0) = eps;
                            postfixshuffPDF1115(postfixshuffPDF1115 == 0) = eps;
                            postfixshuffPDF1620(postfixshuffPDF1620 == 0) = eps;
                            postfixshuffPDF2125(postfixshuffPDF2125 == 0) = eps;
                            postfixshuffPDF(postfixshuffPDF == 0) = eps;
                            
                            postfixshuffPDF15 = postfixshuffPDF15/sum(sum(postfixshuffPDF15));
                            postfixshuffPDF610 = postfixshuffPDF610/sum(sum(postfixshuffPDF610));
                            postfixshuffPDF1115 = postfixshuffPDF1115/sum(sum(postfixshuffPDF1115));
                            postfixshuffPDF1620 = postfixshuffPDF1620/sum(sum(postfixshuffPDF1620));
                            postfixshuffPDF2125 = postfixshuffPDF2125/sum(sum(postfixshuffPDF2125));
                            postfixshuffPDF = postfixshuffPDF/sum(sum(postfixshuffPDF));
                            
                            KLnorm(index,1,t) = sum(sum(log2(prefixPDF15./postfixPDF15).*prefixPDF15))...
                                +sum(sum(log2(postfixPDF15./prefixPDF15).*postfixPDF15));
                            if ii >= 10;
                                KLnorm(index,2,t) = sum(sum(log2(prefixPDF610./postfixPDF610).*prefixPDF610))...
                                    +sum(sum(log2(postfixPDF610./prefixPDF610).*postfixPDF610));
                            end
                            if ii >= 15;
                                KLnorm(index,3,t) = sum(sum(log2(prefixPDF1115./postfixPDF1115).*prefixPDF1115))...
                                    +sum(sum(log2(postfixPDF1115./prefixPDF1115).*postfixPDF1115));
                            end
                            if ii >= 20
                                KLnorm(index,4,t) = sum(sum(log2(prefixPDF1620./postfixPDF1620).*prefixPDF1620))...
                                    +sum(sum(log2(postfixPDF1620./prefixPDF1620).*postfixPDF1620));
                            end
                            if ii >= 25
                                KLnorm(index,5,t) = sum(sum(log2(prefixPDF2125./postfixPDF2125).*prefixPDF2125))...
                                    +sum(sum(log2(postfixPDF2125./prefixPDF2125).*postfixPDF2125));
                            end
                            KLnorm(index,6,t) = sum(sum(log2(prefixPDF./postfixPDF).*prefixPDF))...
                                +sum(sum(log2(postfixPDF./prefixPDF).*postfixPDF));
                            KLshuff(index,1,t) = sum(sum(log2(prefixshuffPDF15./postfixshuffPDF15).*prefixshuffPDF15))...
                                +sum(sum(log2(postfixshuffPDF15./prefixshuffPDF15).*postfixshuffPDF15));
                            if ii >= 10;
                                KLshuff(index,2,t) = sum(sum(log2(prefixshuffPDF610./postfixshuffPDF610).*prefixshuffPDF610))...
                                    +sum(sum(log2(postfixshuffPDF610./prefixshuffPDF610).*postfixshuffPDF610));
                            end
                            if ii >= 15;
                                KLshuff(index,3,t) = sum(sum(log2(prefixshuffPDF1115./postfixshuffPDF1115).*prefixshuffPDF1115))...
                                    +sum(sum(log2(postfixshuffPDF1115./prefixshuffPDF1115).*postfixshuffPDF1115));
                            end
                            if ii >= 20
                                KLshuff(index,4,t) = sum(sum(log2(prefixshuffPDF1620./postfixshuffPDF1620).*prefixshuffPDF1620))...
                                    +sum(sum(log2(postfixshuffPDF1620./prefixshuffPDF1620).*postfixshuffPDF1620));
                            end
                            if ii >= 25
                                KLshuff(index,5,t) = sum(sum(log2(prefixshuffPDF2125./postfixshuffPDF2125).*prefixshuffPDF2125))...
                                    +sum(sum(log2(postfixshuffPDF2125./prefixshuffPDF2125).*postfixshuffPDF2125));
                            end
                            KLshuff(index,6,t) = sum(sum(log2(prefixshuffPDF./postfixshuffPDF).*prefixshuffPDF))...
                                +sum(sum(log2(postfixshuffPDF./prefixshuffPDF).*postfixshuffPDF));
                        end
                    end
                end
            end
        end
    end
end
save('Memory-and-Salience')
%%
load(['C:\Users\skoenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\'...
    'CombinedSalienceStatistics.mat'],'allstatistics','allshuffled');
novstd = zeros(1,25);
repstd = zeros(1,25);
pvalus = zeros(1,25);
zpvalues = zeros(2,25);
shuffled = allshuffled(:,1,1,:);
shuffled(isnan(shuffled)) = [];
for i = 1:25
    novs = novsalfix(:,i,:);
    novs(isnan(novs)) = [];
    novstd(i) = std(novs)/sqrt(length(novs));
    reps = repsalfix(:,i,:);
    reps(isnan(reps)) = [];
    repstd(i) = std(reps)/sqrt(length(reps));
    [~,p] = kstest2(novs,reps);
    pvalues(i) = p;
    [~,p] = ztest(novs,mean(shuffled),std(shuffled));
    zpvalues(1,i) = p;
    [~,p] = ztest(reps,mean(shuffled),std(shuffled));
    zpvalues(2,i) = p;
end

figure
title('Salience at fixations during Novel Viewing and Repeated Viewing')
hold on
errorbar(nanmean(nanmean(novsalfix(:,1:25,:),1),3),novstd);
errorbar(nanmean(nanmean(repsalfix(:,1:25,:),1),3),repstd,'r');
plot(1:25,allstatistics.confidenceintervals(1,1,1)*ones(1,25),'k--');
for i = 1:25;
    if pvalues(i) < 0.05
        plot(i,nanmean(nanmean(novsalfix(:,i,:)))+novstd(i)+0.01,'*r');
    end
end
hold off
legend('Salience @ Novel Fixations','Salience @ Repeated Fixations',...
    'Chance Salience','p < 0.05','location','northeastoutside')
xlim([0,25])


novsal = novsalfix(:,1:25,:);
novsal = novsal(1:end);
repsal = repsalfix(:,1:25,:);
repsal = repsal(1:end);
figure
bar([mean(nanmean(nanmean(novsalfix(:,1:25,:),1),3)),mean(nanmean(nanmean(repsalfix(:,1:25,:),1),3))])
hold on
errorbar([nanmean(novsal) nanmean(repsal)],[nanstd(novsal)/sqrt(sum(~isnan(novsal))) nanstd(novsal)/sqrt(sum(~isnan(novsal)))],'.k')
box off
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Salience')
[~,allsal_pval] = ttest2(novsal,repsal)
%%
clr = ['rgbm'];
Fivepointshuff = KLshuff(:,1:5,:);
Fivepointshuff(isnan(Fivepointshuff)) = [];
TwentyFivepointshuff = KLshuff(:,6,:);
TwentyFivepointshuff(isnan(TwentyFivepointshuff)) = [];
figure
hold all
h(1) = area(ones(1,5)*mean(Fivepointshuff)+std(Fivepointshuff));
set(h(1),'FaceColor',[0.7 0.7 0.7]);
set(h(1),'EdgeColor',[0.7 0.7 0.7]);
h(2) = area(ones(1,5)*mean(Fivepointshuff)-std(Fivepointshuff));
set(h(2),'FaceColor','w');
set(h(2),'EdgeColor','w');
h(3) = area(5:6,ones(1,2)*mean(Fivepointshuff)+std(Fivepointshuff));
set(h(3),'FaceColor',[0.7 0.7 0.7]);
set(h(3),'EdgeColor',[0.7 0.7 0.7]);
h(4) = area(5:6,ones(1,2)*mean(TwentyFivepointshuff)-std(TwentyFivepointshuff));
set(h(4),'FaceColor','w');
set(h(4),'EdgeColor','w');
for t = 1:length(tags);
    p(t) = plot(nanmean(nanmean(KLnorm(:,:,t),3)),clr(t));
end

legend([h(1) p],[{'Chance Distance'} tags],'location','Northeastoutside')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'XTickLabel',{'1-5','6-10','11-15','16-20','21-25','All'})
ylabel('KL Divergence-Distance (Bits)')
ylim([0 70])
%% For Timmy Only
figure
subplot(1,2,1)
title('KL Distances for Timmy')
hold on
plot(nanmean(nanmean(KLnorm(1:144,:,1),3)),'k')
plot(nanmean(nanmean(KLnorm(145:end,:,1),3)),'r')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'XTickLabel',{'1-5','6-10','11-15','16-20','21-25','All'})
ylabel('KL Divergence-Distance (Bits)')
ylim([0 70])
legend('pre-lesion','post-lesion')

subplot(1,2,2)
title('Average Salience at fixations 1-35')
bar([nanmean(nanmean(novsalfix(1:144,1:35,1))),...
    nanmean(nanmean(novsalfix(145:end,1:35,1))),...
    nanmean(nanmean(repsalfix(1:144,1:35,1))),...
    nanmean(nanmean(repsalfix(145:end,1:35,1)))]');
set(gca,'XTick',[1 2 3 4])
set(gca,'XTickLabel',{'pre-Novel','post-Novel','pre-Rep','post-Rep'})
ylabel('Normalized Salience')

novvals1 = novsalfix(1:144,1:35,1);
novvals1(isnan(novvals1)) = [];
novvals2 = novsalfix(145:288,1:35,1);
novvals2(isnan(novvals2)) = [];
[~,novpval] = ttest2(novvals1,novvals2);

repvals1 = repsalfix(1:144,1:35,1);
repvals1(isnan(repvals1)) = [];
repvals2 = repsalfix(145:288,1:35,1);
repvals2(isnan(repvals2)) = [];
[~,reppval] = ttest2(repvals1,repvals2);

[~,prenovrep] = ttest2(novvals1,repvals1);
[~,postnovrep] = ttest2(novvals2,repvals2);
%%
figure
suptitle('Novel viewing Fixation PDFs')
subplot(2,3,1)
imagesc(imfilter(prebiasPDF15,f))
title('Fixations # 1-5')
axis off
subplot(2,3,2)
imagesc(imfilter(prebiasPDF610,f))
title('Fixations # 6-10')
axis off
subplot(2,3,3)
imagesc(imfilter(prebiasPDF1115,f))
title('Fixations # 11-15')
axis off
subplot(2,3,4)
imagesc(imfilter(prebiasPDF1620,f))
title('Fixations # 16-20')
axis off
subplot(2,3,5)
imagesc(imfilter(prebiasPDF2125,f))
title('Fixations # 21- 25')
axis off
subplot(2,3,6)
imagesc(imfilter(prebiasPDF,f))
title('All fixations')
axis off

figure
suptitle('Repeated viewiing Fixation PDFs')
subplot(2,3,1)
imagesc(imfilter(postbiasPDF15,f))
title('Fixations # 1-5')
axis off
subplot(2,3,2)
imagesc(imfilter(postbiasPDF610,f))
title('Fixations # 6-10')
axis off
subplot(2,3,3)
imagesc(imfilter(postbiasPDF1115,f))
title('Fixations # 11-15')
axis off
subplot(2,3,4)
imagesc(imfilter(postbiasPDF1620,f))
title('Fixations # 16-20')
axis off
subplot(2,3,5)
imagesc(imfilter(postbiasPDF2125,f))
title('Fixations # 21- 25')
axis off
subplot(2,3,6)
imagesc(imfilter(postbiasPDF,f))
title('All fixations')
axis off
%%
figure
suptitle('Random "Novel" viewing Fixation PDFs')
subplot(2,3,1)
imagesc(imfilter(prebiasshuffPDF15,f))
title('Fixations # 1-5')
axis off
subplot(2,3,2)
imagesc(imfilter(prebiasshuffPDF610,f))
title('Fixations # 6-10')
axis off
subplot(2,3,3)
imagesc(imfilter(prebiasshuffPDF1115,f))
title('Fixations # 11-15')
axis off
subplot(2,3,4)
imagesc(imfilter(prebiasshuffPDF1620,f))
title('Fixations # 16-20')
axis off
subplot(2,3,5)
imagesc(imfilter(prebiasshuffPDF2125,f))
title('Fixations # 21- 25')
axis off
subplot(2,3,6)
imagesc(imfilter(prebiasshuffPDF,f))
title('All fixations')
axis off

figure
suptitle('Random "Repeated" viewing Fixation PDFs')
subplot(2,3,1)
imagesc(imfilter(postbiasshuffPDF15,f))
title('Fixations # 1-5')
axis off
subplot(2,3,2)
imagesc(imfilter(postbiasshuffPDF610,f))
title('Fixations # 6-10')
axis off
subplot(2,3,3)
imagesc(imfilter(postbiasshuffPDF1115,f))
title('Fixations # 11-15')
axis off
subplot(2,3,4)
imagesc(imfilter(postbiasshuffPDF1620,f))
title('Fixations # 16-20')
axis off
subplot(2,3,5)
imagesc(imfilter(postbiasshuffPDF2125,f))
title('Fixations # 21- 25')
axis off
subplot(2,3,6)
imagesc(imfilter(postbiasshuffPDF,f))
title('All fixations')
axis off
%%
binned_postbias15= bin2(imfilter(postbiasPDF15,f),binsize,binsize);
binned_postbias610= bin2(imfilter(postbiasPDF610,f),binsize,binsize);
binned_postbias1115= bin2(imfilter(postbiasPDF1115,f),binsize,binsize);
binned_postbias1620= bin2(imfilter(postbiasPDF1620,f),binsize,binsize);
binned_postbias2125= bin2(imfilter(postbiasPDF2125,f),binsize,binsize);
binned_postbias= bin2(imfilter(postbiasPDF,f),binsize,binsize);

binned_postbias15(binned_postbias15 == 0) = eps;
binned_postbias15 = binned_postbias15/sum(sum(binned_postbias15));
binned_postbias610(binned_postbias610 == 0) = eps;
binned_postbias610 = binned_postbias610/sum(sum(binned_postbias610));
binned_postbias1115(binned_postbias1115 == 0) = eps;
binned_postbias1115 = binned_postbias1115/sum(sum(binned_postbias1115));
binned_postbias1620(binned_postbias1620 == 0) = eps;
binned_postbias1620 = binned_postbias1620/sum(sum(binned_postbias1620));
binned_postbias2125(binned_postbias2125 == 0) = eps;
binned_postbias2125 = binned_postbias2125/sum(sum(binned_postbias2125));
binned_postbias(binned_postbias == 0) = eps;
binned_postbias = binned_postbias/sum(sum(binned_postbias));

binned_prebias15= bin2(imfilter(prebiasPDF15,f),binsize,binsize);
binned_prebias610= bin2(imfilter(prebiasPDF610,f),binsize,binsize);
binned_prebias1115= bin2(imfilter(prebiasPDF1115,f),binsize,binsize);
binned_prebias1620= bin2(imfilter(prebiasPDF1620,f),binsize,binsize);
binned_prebias2125= bin2(imfilter(prebiasPDF2125,f),binsize,binsize);
binned_prebias= bin2(imfilter(prebiasPDF,f),binsize,binsize);

binned_prebias15(binned_prebias15 == 0) = eps;
binned_prebias15 = binned_prebias15/sum(sum(binned_prebias15));
binned_prebias610(binned_prebias610 == 0) = eps;
binned_prebias610 = binned_prebias610/sum(sum(binned_prebias610));
binned_prebias1115(binned_prebias1115 == 0) = eps;
binned_prebias1115 = binned_prebias1115/sum(sum(binned_prebias1115));
binned_prebias1620(binned_prebias1620 == 0) = eps;
binned_prebias1620 = binned_prebias1620/sum(sum(binned_prebias1620));
binned_prebias2125(binned_prebias2125 == 0) = eps;
binned_prebias2125 = binned_prebias2125/sum(sum(binned_prebias2125));
binned_prebias(binned_prebias == 0) = eps;
binned_prebias = binned_prebias/sum(sum(binned_prebias));

avgKL(1) = sum(sum(log2(binned_prebias15./binned_postbias15).*binned_prebias15))...
    +sum(sum(log2(binned_postbias15./binned_prebias15).*binned_postbias15));
avgKL(2) = sum(sum(log2(binned_prebias610./binned_postbias610).*binned_prebias610))...
    +sum(sum(log2(binned_postbias610./binned_prebias610).*binned_postbias610));
avgKL(3) = sum(sum(log2(binned_prebias1115./binned_postbias1115).*binned_prebias1115))...
    +sum(sum(log2(binned_postbias1115./binned_prebias1115).*binned_postbias1115));
avgKL(4) = sum(sum(log2(binned_prebias1620./binned_postbias1620).*binned_prebias1620))...
    +sum(sum(log2(binned_postbias1620./binned_prebias1620).*binned_postbias1620));
avgKL(5) = sum(sum(log2(binned_prebias2125./binned_postbias2125).*binned_prebias2125))...
    +sum(sum(log2(binned_postbias2125./binned_prebias2125).*binned_postbias2125));
avgKL(6) = sum(sum(log2(binned_prebias ./binned_postbias ).*binned_prebias ))...
    +sum(sum(log2(binned_postbias ./binned_prebias ).*binned_postbias ));

figure
plot(avgKL)
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'XTickLabel',{'1-5','6-10','11-15','16-20','21-25','All'})
ylabel('KL Divergence-Distance (Bits)')
title('Distance between the average Fixation PDFs')
%%
binned_postbiasshuff15= bin2(imfilter(postbiasshuffPDF15,f),binsize,binsize);
binned_postbiasshuff610= bin2(imfilter(postbiasshuffPDF610,f),binsize,binsize);
binned_postbiasshuff1115= bin2(imfilter(postbiasshuffPDF1115,f),binsize,binsize);
binned_postbiasshuff1620= bin2(imfilter(postbiasshuffPDF1620,f),binsize,binsize);
binned_postbiasshuff2125= bin2(imfilter(postbiasshuffPDF2125,f),binsize,binsize);
binned_postbiasshuff= bin2(imfilter(postbiasshuffPDF,f),binsize,binsize);

binned_postbiasshuff15(binned_postbiasshuff15 == 0) = eps;
binned_postbiasshuff15 = binned_postbiasshuff15/sum(sum(binned_postbiasshuff15));
binned_postbiasshuff610(binned_postbiasshuff610 == 0) = eps;
binned_postbiasshuff610 = binned_postbiasshuff610/sum(sum(binned_postbiasshuff610));
binned_postbiasshuff1115(binned_postbiasshuff1115 == 0) = eps;
binned_postbiasshuff1115 = binned_postbiasshuff1115/sum(sum(binned_postbiasshuff1115));
binned_postbiasshuff1620(binned_postbiasshuff1620 == 0) = eps;
binned_postbiasshuff1620 = binned_postbiasshuff1620/sum(sum(binned_postbiasshuff1620));
binned_postbiasshuff2125(binned_postbiasshuff2125 == 0) = eps;
binned_postbiasshuff2125 = binned_postbiasshuff2125/sum(sum(binned_postbiasshuff2125));
binned_postbiasshuff(binned_postbiasshuff == 0) = eps;
binned_postbiasshuff = binned_postbiasshuff/sum(sum(binned_postbiasshuff));

binned_prebiasshuff15= bin2(imfilter(prebiasshuffPDF15,f),binsize,binsize);
binned_prebiasshuff610= bin2(imfilter(prebiasshuffPDF610,f),binsize,binsize);
binned_prebiasshuff1115= bin2(imfilter(prebiasshuffPDF1115,f),binsize,binsize);
binned_prebiasshuff1620= bin2(imfilter(prebiasshuffPDF1620,f),binsize,binsize);
binned_prebiasshuff2125= bin2(imfilter(prebiasshuffPDF2125,f),binsize,binsize);
binned_prebiasshuff= bin2(imfilter(prebiasshuffPDF,f),binsize,binsize);

binned_prebiasshuff15(binned_prebiasshuff15 == 0) = eps;
binned_prebiasshuff15 = binned_prebiasshuff15/sum(sum(binned_prebiasshuff15));
binned_prebiasshuff610(binned_prebiasshuff610 == 0) = eps;
binned_prebiasshuff610 = binned_prebiasshuff610/sum(sum(binned_prebiasshuff610));
binned_prebiasshuff1115(binned_prebiasshuff1115 == 0) = eps;
binned_prebiasshuff1115 = binned_prebiasshuff1115/sum(sum(binned_prebiasshuff1115));
binned_prebiasshuff1620(binned_prebiasshuff1620 == 0) = eps;
binned_prebiasshuff1620 = binned_prebiasshuff1620/sum(sum(binned_prebiasshuff1620));
binned_prebiasshuff2125(binned_prebiasshuff2125 == 0) = eps;
binned_prebiasshuff2125 = binned_prebiasshuff2125/sum(sum(binned_prebiasshuff2125));
binned_prebiasshuff(binned_prebiasshuff == 0) = eps;
binned_prebiasshuff = binned_prebiasshuff/sum(sum(binned_prebiasshuff));

avgKLshuff(1) = sum(sum(log2(binned_prebiasshuff15./binned_postbiasshuff15).*binned_prebiasshuff15))...
    +sum(sum(log2(binned_postbiasshuff15./binned_prebiasshuff15).*binned_postbiasshuff15));
avgKLshuff(2) = sum(sum(log2(binned_prebiasshuff610./binned_postbiasshuff610).*binned_prebiasshuff610))...
    +sum(sum(log2(binned_postbiasshuff610./binned_prebiasshuff610).*binned_postbiasshuff610));
avgKLshuff(3) = sum(sum(log2(binned_prebiasshuff1115./binned_postbiasshuff1115).*binned_prebiasshuff1115))...
    +sum(sum(log2(binned_postbiasshuff1115./binned_prebiasshuff1115).*binned_postbiasshuff1115));
avgKLshuff(4) = sum(sum(log2(binned_prebiasshuff1620./binned_postbiasshuff1620).*binned_prebiasshuff1620))...
    +sum(sum(log2(binned_postbiasshuff1620./binned_prebiasshuff1620).*binned_postbiasshuff1620));
avgKLshuff(5) = sum(sum(log2(binned_prebiasshuff2125./binned_postbiasshuff2125).*binned_prebiasshuff2125))...
    +sum(sum(log2(binned_postbiasshuff2125./binned_prebiasshuff2125).*binned_postbiasshuff2125));
avgKLshuff(6) = sum(sum(log2(binned_prebiasshuff ./binned_postbiasshuff ).*binned_prebiasshuff ))...
    +sum(sum(log2(binned_postbiasshuff ./binned_prebiasshuff ).*binned_postbiasshuff ));

figure
plot(avgKLshuff)
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'XTickLabel',{'1-5','6-10','11-15','16-20','21-25','All'})
ylabel('KL Divergence-Distance (Bits)')
title('Distance between the average Shuffled Fixation PDFs')
%%
pvalues = ones(1,6);
means = ones(2,6);
stds = ones(2,6);
for i = 1:6;
    shuff = KLshuff(:,i,:);
    norm = KLnorm(:,i,:);
    shuff(isnan(shuff)) = [];
    norm(isnan(norm)) = [];
    [~,p] = ttest2(norm,shuff);
    pvalues(i) = p;
    means(:,i) = [mean(norm); mean(shuff)];
    stds(:,i) = [std(norm); std(shuff)];
end
    