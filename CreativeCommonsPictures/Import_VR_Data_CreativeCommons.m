data_dir = '\\research.wanprc.org\research\Buffalo Lab\Virtual Navigation\UnityVR\EyeBackup\';
img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\CreativeCommonsPictures\Renamed\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\CreativeCommonsPictures\Renamed\Figures\';
eyedata_dir = img_dir;

% data_file = '16_06_21_14_32';
% 
% eye_cal_file = ['eye_cal2_' data_file]; %contains time, x, and y data
% time_cal_file =['time_cal2_' data_file]; %contains event and tieyeing data
% 
% 
% %---Read in Eye and Event Data--%
% [eye1 eye2 eye3] = textread([data_dir eye_cal_file],'%f%f%f','headerlines',1,'delimiter',',','whitespace','\n');
% 
% try %sometimes eye data file gets cutoff at end?
%     eye = [eye1 eye2 eye3];
% catch
%     eye = [eye1(1:end-1) eye2 eye3];
% end
% 
% t0 = eye(1,1); %time 0
% eye(:,1) = eye(:,1)-t0; %make time relative to 0 otherwise relative to 1970
% 
% time = NaN(1,1000); %prealocate space for time of events
% square_pos = NaN(1000,2); %prealocate space for square_position of calibration squares
% events = NaN(1,1000); %prealocate space for events
% % 1: square square_position, 2: square on, 3:reward, 4:break fixation
% % 5: clrchng, 6: square off
% 
% fid = fopen([data_dir time_cal_file]);
% linecount = 1;
% tline = fgets(fid); %skip line 1: timestaeyep, events
% tline = fgets(fid); %skip line 2: [timestaeyep1] start collecting eye data
% image_name = {};
% img_count = 1;
% while ischar(tline)
%     tline = fgets(fid); %get next line
%     if tline ~= -1; %end of file is noted by -1
%         C = textscan(tline,'%s'); %parse line by spaces and comas
%         time(linecount) = str2double(C{1}{1});
%         if strcmpi(C{1}{2},'Square')
%             if strcmpi(C{1}{3},'Position,') %gives position of square in pixels
%                 events(linecount) = 1;
%                 square_pos(linecount,1) = str2double(C{1}{4});
%                 square_pos(linecount,2) = str2double(C{1}{5});
%             elseif strcmpi(C{1}{3},'moved') || strcmpi(C{1}{3},'on') %square moved to different location and displayed
%                 events(linecount) = 2;
%             elseif strcmpi(C{1}{3},'dims') %color change occured
%                 events(linecount) = 5;
%             elseif strcmpi(C{1}{3},'off') %suqre diseappears
%                 events(linecount) = 6;
%             else
%                 disp('unknown events parameter in column 3 of event array') %so we known we didn't miss anything
%             end
%         elseif strcmpi(C{1}{2},'reward') %monkey received reward
%             events(linecount) = 3;
%         elseif strcmpi(C{1}{2},'no') %trial aborted due to break fixation, etc.
%             events(linecount) = 4;
%         elseif ~isempty(strfind(lower(C{1}{2}),'photo')) %image trials
%             if length(C{1}) == 3
%                 if  strcmpi(C{1}{3},'on') %image on
%                     events(linecount) = 7;
%                 elseif strcmpi(C{1}{3},'off') %image off
%                     events(linecount) = 8;
%                 end
%             else
%                 if~isempty(strfind(C{1}{2},'/')) %image name
%                     events(linecount) = 9;
%                     slash = strfind(C{1}{2},'/');
%                     image_name{img_count} = C{1}{2}(slash+1:end);
%                     img_count = img_count+1;
%                 else
%                     disp('Unknown parameter for column 2 for picture trials')  %so we known we didn't miss anything
%                 end
%             end
%         end
%         linecount = linecount+1; %to go to the next line
%     end
% end
% fclose(fid);
% 
% %---Get Set number for images...for later use---%
% if strcmpi(image_name{1}(1),'S')
%     setnum = str2double(image_name{1}(2:3)); %get the set number, assumes set number is the same for all images
%     if setnum < 10
%         img_dir = ['C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\SCM0' num2str(setnum) '\'];
%     else
%         img_dir = ['C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\SCM' num2str(setnum) '\'];
%         
%     end
% else %VR pilot sets
%     setnum = str2double(image_name{1}(4)); %get the set number, assumes set number is the same for all images
%     img_dir = ['C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\VRSCM Images\VRset' num2str(setnum) '\'];
% end
% 
% % may be incorrect
% set_dir = [figure_dir 'VRSet' num2str(setnum) '\'];
% mkdir(set_dir);  %make directory if it doesn't already
% %exist for individual sets for organization purpose
% 
% %---Get the calibration Function---%
% %doing for individual sessions of calibration. Currently assumes stable
% %cailbration within short session. Does not track drift
% if linecount < 1000 %if you didn't use the whole prealocated space remove extra rows
%     time = time(1:linecount-1);
%     events = events(1:linecount-1);
%     square_pos = square_pos(1:linecount-1,:);
% end
% time = time-t0; %use t0 to ensure time is the saeyee but should be
% 
% % get unique calibration points
% xsquare_poses = unique(square_pos(:,1)); xsquare_poses(isnan(xsquare_poses)) = [];
% ysquare_poses = unique(square_pos(:,2)); ysquare_poses(isnan(ysquare_poses)) = [];
% 
% 
% rewards = find(events == 3); %find where a reward occured
% caldat_x = cell(length(xsquare_poses),length(ysquare_poses));
% caldat_y = cell(length(xsquare_poses),length(ysquare_poses));
% square_dims = find(events == 5);
% square_pos_events = find(events == 1);
% for rew = 1:length(rewards);
%     rewardind = rewards(rew); %index of event array in which reward was delivered
%     
%     %trialstarttieye = time(rewardind-3); %square on
%     square_change = square_dims(square_dims < rewardind); %find the last time the color change
%     square_change = square_change(end);
%     clrchngtieye = time(square_change); %look at time in which square changes color (dims)
%     start_data_collection = clrchngtieye - 0.5; %500 eyes
%     
%     %find indices associated with times between 500 ms before color change
%     %and color change
%     eyeind = find(eye(:,1) > start_data_collection & eye(:,1) <= clrchngtieye);
%     
%     %determine x and y eye position for those 500 ms
%     xsquare_pos = eye(eyeind,2);
%     ysquare_pos = eye(eyeind,3);
%     
%     %position of square
%     sq_pos = square_pos_events(square_pos_events < rewardind);
%     sq_pos = sq_pos(end);
%     square_x = square_pos(sq_pos,1);
%     square_y = square_pos(sq_pos,2);
%     
%     xind = find(xsquare_poses == square_x);
%     yind = find(ysquare_poses == square_y);
%     
%     caldat_x{xind,yind} = [caldat_x{xind,yind} mean(xsquare_pos)];
%     caldat_y{xind,yind} = [caldat_y{xind,yind} mean(ysquare_pos)];
% end
% 
% % control are unique combinations x and y position of square in pixel coordinates
% control_x = [];
% control_y = [];
% for xi = 1:length(xsquare_poses);
%     for yi = 1:length(ysquare_poses);
%         control_x = [control_x; 1.5*xsquare_poses(xi)];
%         control_y = [control_y; ysquare_poses(yi)];
%     end
% end
% 
% % input are estimated position of the color change squares in mV (the ouptu
% % of eyescan).
% input_x = [];
% input_y = [];
% for xi = 1:length(xsquare_poses);
%     for yi = 1:length(ysquare_poses);
%         input_x = [input_x; nanmean(caldat_x{xi,yi})];
%         input_y = [input_y; nanmean(caldat_y{xi,yi})];
%     end
% end
% 
% %tform automatically transforms estimated eye positions in mV into pixel
% %coordinates in  a "magical" way.
% tform = cp2tform([control_x control_y], [input_x input_y],'polynomial',4);
% tform.forward_fcn = tform.inverse_fcn;
% 
% 
% %---Create Plot of Calibration Quality---%
% %Visualize position of estimated square position (blue *) compared to
% %actual position of color change squares (red +)
% [xp,yp] = tformfwd(tform,input_x,input_y);%transform inputs from mV to pixels to plot
% figure
% hold on
% for xi = 1:length(xsquare_poses);
%     for yi = 1:length(ysquare_poses);
%         plot(1.5*xsquare_poses(xi),ysquare_poses(yi),'r+')
%     end
% end
% for xi = 1:length(xp);
%     plot(xp(xi),yp(xi),'b*')
% end
% hold off
% save_and_close_fig(set_dir,['Calibration for ' data_file])
% 
% 
% %---Recalibrate Collected Eye Data---%
% eyedata = {};
% img_on = find(events == 7);
% img_off = find(events == 8);
% for img = 1:length(image_name);
%     eye_ind = find(eye(:,1) > time(img_on(img)) & eye(:,1) <= time(img_off(img)));
%     
%     [xp,yp] = tformfwd(tform,eye(eye_ind,2)',eye(eye_ind,3)');%transform inputs from mV to pixels to plot
%     eyedat{img} = [xp;yp];
% end
% 
% %---Get image numbers and the order they appeared in---%
% imgnum = NaN(1,108);
% for img = 1:length(image_name);
%     dot = strfind(image_name{img},'.');
%     if double(image_name{img}(dot-1)) > 64 %character str
%          imgnum(img) = str2double(image_name{img}(1:dot-2));
%     else
%         imgnum(img) = str2double(image_name{img}(1:dot-1));
%     end
% end
% 
% imgnum(imgnum == 0) = NaN;%don't know why zero shows up
% 
% pairings = NaN(2,54);
% for img = min(imgnum):max(imgnum) %assumes images don't cross over sets
%     ind = find(imgnum == img);
%     if length(ind) == 2 %otherwise don't care
%         pairings(1,img) = ind(1); %novel
%         pairings(2,img) = ind(2); %repeat
%     end
% end
% 
% for eye = 1:length(eyedat)
%     x = eyedat{eye}(1,:);
%     y = eyedat{eye}(2,:);
%     x = x+400;
%     y = y+300;
%     
%     y(x < -50) = [];
%     x(x < -50) = [];
%     x(y < -50) = [];
%     y(y < -50) = [];
%     
%     y(x > 850) = [];
%     x(x > 850) = [];
%     x(y > 650) = [];
%     y(y > 650) = [];
%     
%     eyedat{eye} = [x;y];
% end
% 
%     fixationstats = ClusterFixation_Short(eyedat);
% save([eyedata_dir data_file '-fixation.mat'],'image_name','pairings',...
%     'imgnum','setnum','eyedat','fixationstats')
%%
% tags = {'MP','TT','JN','IW'};
% f = fspecial('gaussian',[256,256],24);
% imageX = 800; imageY = 600;
% %---Plot Eye Data on Corresponding Figure---%
% for img = 26%1:size(pairings,2);
%     im1 = imread([img_dir image_name{pairings(1,img)}]);
%     
%     load([img_dir image_name{pairings(1,img)}(1:end-4) '-saliencemap.mat'],'fullmap');
%     
%     allBCRW = zeros(imageY,imageX);
%     for t = 1:length(tags)
%         load([img_dir 'BCRW IOR TAU 17\' tags{t} '-' num2str(img) '-BCRW.mat'],'fixations')
%         allBCRW = allBCRW+fixations;
%     end
%     allBCRW = imfilter(allBCRW,f);
%     
%     x = eyedat{pairings(1,img)}(1,:);
%     y = 600-eyedat{pairings(1,img)}(2,:);
%     
%     
%     figure
%     subplot(2,2,1)
%     imshow(im1)
%     hold on
%     plot(x,y,'w')
%     hold off
%     
%     subplot(2,2,2)
%     imagesc(fullmap)
%     hold on
%     plot(x,y,'w')
%     hold off
%     axis off
%     axis equal
%     colormap('jet')
%     
%     subplot(2,2,3)
%     imagesc(allBCRW)
%     hold on
%     plot(x,y,'w')
%     hold off
%     axis off
%     axis equal
%     colormap('jet')
%     
%     save_and_close_fig(figure_dir,[num2str(img) '_maps_plus_eyetrace'])
% end
%%
imgind = 36; %for image 26 novel

allBCRW = zeros(imageY,imageX);
for t = 1:length(tags)
    load([img_dir 'BCRW IOR TAU 17\' tags{t} '-26-BCRW.mat'],'fixations')
    allBCRW = allBCRW+fixations;
end
f = fspecial('gaussian',[256,256],24);
allBCRW = imfilter(allBCRW,f);
load([img_dir '26-saliencemap.mat'],'fullmap');

fixations = round(fixationstats{36}.fixations);


BCRW_vals = [];
rBCRW_vals = [];
sal_vals = []; 
rsal_vals = [];
for f = 1:length(fixations)
    
    if fixations(1,f) < 1 || fixations(1,f) > imageX || fixations(2,f) < 1 || fixations(2,f) > imageY
       continue 
    end
    
    BCRW_vals = [BCRW_vals allBCRW(imageY-fixations(2,f),fixations(1,f))]; 
    sal_vals =  [sal_vals allBCRW(imageY-fixations(2,f),fixations(1,f))]; 
    
    ry = randi(600); ry(ry < 1) = 1;
    rx = randi(800); rx(rx < 1) = 1;
    rsal_vals = [rsal_vals fullmap(ry,rx)];
    rBCRW_vals = [rBCRW_vals allBCRW(ry,rx)];
end

thresh = 0:0.01:1;
TP = NaN(2,length(thresh)); %True positive
FA = NaN(2,length(thresh)); %False alarm
len = length(sal_vals);
for ii = 1:length(thresh)
    TP(1,ii) = sum(sal_vals > thresh(ii))/len;
    FA(1,ii) = sum(rsal_vals > thresh(ii))/len;
    TP(2,ii) = sum(BCRW_vals > thresh(ii))/len;
    FA(2,ii) = sum(rBCRW_vals > thresh(ii))/len;
end

figure
plot(FA(1,:),TP(1,:),'b')
hold on
plot(FA(2,:),TP(2,:),'g')
plot([0 1],[0 1],'k--')
hold off
legend('Salience','BCRW')

-trapz(FA(1,:),TP(1,:))
-trapz(FA(2,:),TP(2,:))



