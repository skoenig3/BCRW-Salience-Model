image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800;
imageY = 600;

IOR_dir = 'BCRW IOR TAU ';
IOR_taus = [0 50 35 25 17 12 7 3 1];


distance_btwn_fx = cell(length(IOR_taus),length(image_sets));
all_dist = cell(1,length(IOR_taus));
for imset = 1:length(image_sets);
    for IOR = 1:length(IOR_taus);
        distance_btwn_fix{IOR,imset} = NaN(36*4,75);
        count = 1;
        for i = 1:36;
            for t = 1:length(tags)
                load([image_dir '\' image_sets{imset} '\' IOR_dir num2str(IOR_taus(IOR)) '\' tags{t} '-' num2str(i) '-BCRW.mat'])
                disp(['Image set-' image_sets{imset} ' IOR_tau = ' num2str(IOR_taus(IOR)) ...
                    ' Image# ' num2str(i) ' Monkey ' tags{t}])
                temp_dist = NaN(size(fixationtimes,1),75);
                for ii = 1:size(fixationtimes,1);
                    img = zeros(imageY,imageX);
                    tind = find(fixationtimes(ii,:,1) > 0);
                    for iii = 2:length(tind)-1
                        x1 = fixationtimes(ii,tind(iii),1);
                        y1 = fixationtimes(ii,tind(iii),2);
                        
                        x2 = fixationtimes(ii,tind(iii+1),1);
                        y2 = fixationtimes(ii,tind(iii+1),2);
                        temp_dist(ii,iii) = sqrt((x1-x2).^2+(y1-y2).^2);
                        
%                         if  temp_dist(ii,iii) < 48
%                             disp('now')
%                         end
                    end
                end
                distance_btwn_fix{IOR,imset}(count,:) = nanmean(temp_dist);
                all_dist{IOR} = [all_dist{IOR}; temp_dist(~isnan(temp_dist))];
                count = count+1;
            end
        end
    end
end
%%
figure
hold all
for IOR = 1:length(IOR_taus)
   [n,x] = hist(all_dist{IOR}/24,250);
   plot(x,n/sum(n))
end
legend({'inf','1/50','1/35','1/25','1/17','1/12','1/7','1/3','1'})
   xlim([0 25])
grid minor
xlabel('Distance between fixations (dva)')
ylabel('Proportion')
%%
prop = [];
for IOR = 1:length(IOR_taus)
    d = all_dist{IOR};
    d = d(1:end)/24;
    prop(IOR) = sum(d < 5)/sum(~isnan(d));
end
prop