function run_BCRWCF_param_sweep(CombinedViewingBehaviorData,saliencemapfile,...
    tagname,imageX,imageY,plotoptions,IOR_area,IOR_tau,border_buffer,border_sacdist,...
    imagesetfolder)
% Function is essentially the same as run_BCRCF.m except tweaked to handle
% a parameter sweep. See run_BCRWCF.m for more details.

warning off all
if nargin < 11
    error(['Not enough inputs: function requires CombinedViewingBehaviorData,'...
        'tagname saliencemapfile, imageX, imageY, and plot options' ....
        'and IOR_area, IOR_tau, border_buffer, border_sacdist,imagesetfolder'])
end

data = CombinedViewingBehaviorData;
distance = data.distanceprofile;
probdst = histc(distance,1:max(max(distance)),1);
distCDF = cumsum(probdst,1);
sCDF = sum(probdst,1);
distCDF = bsxfun(@rdivide,distCDF,sCDF);
sacend = data.mediansac;
sacdistance = distCDF(:,1:sacend);
fixdistance = distCDF(:,sacend+1:end);
persistence.sac = nanmean(data.persistence(:,1:sacend));
persistence.fix = nanmean(data.persistence(:,sacend+1:end));
sacduration = data.sacduration(1:end); sacduration(isnan(sacduration)) = [];
probsacduration = hist(sacduration,max(sacduration));
sacdurationCDF = cumsum(probsacduration)/sum(probsacduration);
fixduration = data.fixduration(1:end); fixduration(isnan(fixduration)) = [];
probfixduration = hist(fixduration,max(fixduration));
fixdurationCDF = cumsum(probfixduration)/sum(probfixduration);
sacangle = data.sacangle; sacangle = sacangle(1:end); sacangle(isnan(sacangle))= [];
sacangle = sacangle*180/pi+180;
nang = (0:360);
probsacangle = hist(sacangle,nang);
sacangleCDF = cumsum(probsacangle)/sum(probsacangle);
load(saliencemapfile,'fullmap');
saliencemap = fullmap;
dash = strfind(saliencemapfile,'-');
dash = dash(1);
clear CombinedViewingBehaviorData data fullmap distCDF

filt = fspecial('gauss',128,12);
% saliencemap(saliencemap < 0.1) = 0; %ignore lowest 10% of salience
saliencemap = imfilter(saliencemap,filt);
[fx fy] = gradient(saliencemap);

nn = 10;%100
trialtime = 10; %10 second simulation
[rr cc] = meshgrid(1:imageX,1:imageY);
dt = 0.005; %5 ms
fixations = zeros(imageY,imageX);
fixationtimes = zeros(nn,trialtime/dt,2);
fixationorder = zeros(600,800);
alltrials = zeros(600,800);
for n = 1:nn;
    if strcmpi(plotoptions.runs,'all')
        figure
        if strcmp(plotoptions.type,'image');
            imagesc(imread([num2str(saliencemapfile(1:dash-1)) '.bmp']))
        else
            imagesc(saliencemap)
        end
        hold on
    end
    fixcount = 1;
    x = imageX/2;
    y = imageY/2;
    xxyy = [[x;y] zeros(2,9)];
    tmr = 0;
    t = 0;
    angold = [];
    
    %select a random saccade duration
    sacdur = find(rand <= sacdurationCDF);
    sacdur = sacdur(1);
    timewarp = round(linspace(1,sacend,sacdur));
    %saccade amplitude and over time is controlled
    sacdist = sacdistance(:,timewarp);
    persac = persistence.sac(timewarp);
    %select a random fixation duration
    fixdur = find(rand <= fixdurationCDF);
    fixdur = fixdur(1);
    timewarp = round(linspace(1,size(fixdistance,2),fixdur));
    %fixation "amplitude" and persistence are controleld over time as well
    fixdist = fixdistance(:,timewarp);
    perfix = persistence.fix(timewarp);
    
    if IOR_tau > 0;
        previous_fixations = NaN(2,1/IOR_tau);
        previous_fixations(:,end) =[x;y];
    end
    fxx = fx;
    fyy = fy;
    C = sqrt((rr-x).^2+(cc-y).^2)<=IOR_area;
    Cind = find(C);
    fxx(Cind) = 0;
    fyy(Cind) = 0;
    while t < trialtime
        if round(tmr*1/dt)+1 <= sacdur;
            dhr = find(sacdist(:,round(tmr*1/dt)+1) >= rand);
            dh = dhr(1);
            b = persac(round(tmr*1/dt)+1)/2;
        else
            dhr = find(fixdist(:,round(tmr*1/dt)-sacdur+1) >= rand);
            dh = dhr(1);
            b = perfix(round(tmr*1/dt)-sacdur+1);
        end
        if round(tmr*1/dt)+1 == sacdur + fixdur
            xy = ceil(mean(xxyy(:,1:5),2));
            if strcmpi(plotoptions.runs,'all')
                plot(xy(1),xy(2),'*k','markersize',6)
            end
            fixations(xy(2),xy(1)) =  fixations(xy(2),xy(1)) + 1;
            fixationorder(xy(2),xy(1)) = fixationorder(xy(2),xy(1)) + fixcount;
            fixationtimes(n,round(t/dt),:) = xy;
            if IOR_tau > 0;
                xyp = previous_fixations(:,1);
                C = sqrt((rr-xyp(1)).^2+(cc-xyp(2)).^2)<=IOR_area;
                fxx(C) = fx(C);
                fyy(C) = fy(C);
                previous_fixations = [previous_fixations(:,2:end) xy];
            end
            C = sqrt((rr-xy(1)).^2+(cc-xy(2)).^2)<=IOR_area;
            Cind = find(C);
            fxx(Cind) = 0;
            fyy(Cind) = 0;
            tmr = 0;
            xxyy = [[x;y] zeros(2,9)];
            angold = [];
            fixcount = fixcount + 1;
            
            sacdur = find(rand <= sacdurationCDF);
            sacdur = sacdur(1);
            timewarp = round(linspace(1,sacend,sacdur));
            sacdist = sacdistance(:,timewarp);
            persac = persistence.sac(timewarp);
            fixdur = find(rand <= fixdurationCDF);
            fixdur = fixdur(1);
            timewarp = round(linspace(1,size(fixdistance,2),fixdur));
            fixdist = fixdistance(:,timewarp);
            perfix = persistence.fix(timewarp);
        end
        if tmr == 0; %just starting this simuliation
            angr = find(sacangleCDF >= rand);
            ang = nang(angr(1));
        else
            if x > imageX-border_buffer || x < border_buffer || ...
                    y > imageY-border_buffer || y < border_buffer %was 10
                if round(tmr*1/dt)+1 <= sacdur
                    [dh ang] = border1(dh,x,y,imageX,imageY,dt,tmr,sacdur,...
                        border_sacdist,border_buffer);
                end
            else
                if abs(fyy(y,x)) == 0 && abs(fxx(y,x)) == 0
                    if round(tmr*1/dt)+1 <= sacdur
                        angh = angold;
                    else
                        angr = randi(length(nang));
                        angh = nang(angr);
                    end
                elseif  abs(fyy(y,x)) > 0 && abs(fxx(y,x)) > 0
                    angh = atand(fyy(y,x)/fxx(y,x));
                    if fxx(y,x) < 0
                        angh = angh + 180;
                    end
                elseif abs(fyy(y,x)) > 0
                    if fyy(y,x) < 0;
                        angh = 90;
                    else
                        angh = 270;
                    end
                elseif abs(fxx(y,x)) > 0
                    if fxx(y,x) < 0;
                        angh = 180;
                    else
                        angh = 0;
                    end
                end
                if angold > 360
                    angold = angold - 360;
                elseif angold < -360
                    angold = angold + 360;
                end
                if abs(angh-angold) > 180
                    if angh < angold
                        angh = angh+360;
                    elseif angold < angh
                        angold = angold+360;
                    end
                end
                ang = angold*(1-b) + b*angh;
            end
        end
        xn = round(x + dh*cos(ang*pi/180));
        yn = round(y + dh*sin(ang*pi/180));
        angold = ang;
        if (xn > imageX || xn < 1 || yn < 1 || yn > imageY)
            [xn yn angold] = border2(x,xn,y,yn,dh,imageX,imageY,dt,tmr,sacdur,...
                border_sacdist,border_buffer);
        end
        if strcmpi(plotoptions.runs,'all')
            plot([x xn],[y yn],'m')
            plot(xn,yn,'.m','markersize',3)
            %             pause(0.01)
        end
        alltrials(yn,xn) = alltrials(yn,xn) + 1;
        x = xn;
        y = yn;
        tmr = tmr + dt;
        t = t+dt;
        xxyy =  [[x;y]  xxyy(:,1:9)];
    end
    if strcmpi(plotoptions.runs,'all')
        close
    end
end

filt = fspecial('gauss',32,6);
if strcmpi(plotoptions.probdens,'all')
    figure
    f = 5*alltrials;
    f = imfilter(f,filt,'replicate');
    f(f<1) = 1;
    imagesc(log(f))
    imagesc(alltrials)
    title('PDF: All Positions')
    figure
    f = 100*fixations;
    f = imfilter(f,filt,'replicate');
    f(f<1) = 1;
    imagesc(log(f))
    title('PDF: Fixations')
    figure,imagesc(saliencemap)
    title('Saliencemap')
end
    function   [dh ang] = border1(dh,x,y,imageX,imageY,dt,tmr,sacdur,...
            border_sacdist,border_buffer)
        if x >= imageX-border_buffer && y >= imageY-border_buffer
            ang = 225;
            dh(dh < border_sacdist*sqrt(2)) = border_sacdist*sqrt(2);
        elseif x >= imageX-border_buffer && y <= border_buffer
            dh(dh < border_sacdist*sqrt(2)) = border_sacdist*sqrt(2);
            ang = 135;
        elseif x <= border_buffer && y >= imageY-border_buffer
            ang = 315;
            dh(dh < border_sacdist*sqrt(2)) = border_sacdist*sqrt(2);
        elseif x <= border_buffer && y <= border_buffer
            ang = 45;
            dh(dh < border_sacdist*sqrt(2)) = border_sacdist*sqrt(2);
        elseif x <= border_buffer
            ang = 0;
            dh(dh < border_sacdist) = border_sacdist;
        elseif y <= border_buffer
            ang = 90;
            dh(dh <  border_sacdist) = border_sacdist;
        elseif x >= imageX-border_buffer
            ang = 180;
            dh(dh < border_sacdist) = border_sacdist;
        elseif y >= imageY-border_buffer
            ang = 270;
            dh(dh < border_sacdist) = border_sacdist;
        end
    end
    function [xn yn angold] = border2(x,xn,y,yn,dh,imageX,imageY,dt,tmr,sacdur,...
            border_sacdis0t,border_buffer)
        xx = [x xn]; yy = [y yn];
        p = polyfit([x xn],[y yn],1);
        if xn > imageX
            y = round(p(1)*imageX + p(2));
            x = imageX;
        elseif xn < 1
            y = round(1*p(1)+p(2));
            x = 1;
        elseif yn < 1
            if any(isinf(p));
                x = x; %#ok
                y = 1;
            else
                x = round((0-p(2))/p(1));
                y = 1;
            end
        elseif yn > imageY
            if  any(isinf(p));
                x = x; %#ok
                y = imageY;
            else
                x = round((imageY-p(2))/p(1));
                y = imageY;
            end
        else
            x = xn; y = yn;
        end
        x(x > imageX) = imageX;
        y(y > imageY) = imageY;
        [dhn angn] = border1(dh,x,y,imageX,imageY,dt,tmr,sacdur,...
            border_sacdist,border_buffer);
        xn = round(x + dhn*cos(angn*pi/180));
        yn = round(y + dhn*sin(angn*pi/180));
        xn(xn > imageX) = imageX; xn(xn < 1) = 1;
        yn(yn > imageY) = imageY; yn(yn < 1) = 1;
        angold = angn;
    end

save([imagesetfolder '\' tagname '-' saliencemapfile(1:dash-1) '-BCRW'],...
    'alltrials', 'fixations','fixationorder','fixationtimes')
%     '-IA_' num2str(IOR_area) '-IT_' num2str(100*IOR_tau) 'BB_' num2str(border_buffer)...
%     '-BS_' num2str(border_sacdist) '.mat'],

end