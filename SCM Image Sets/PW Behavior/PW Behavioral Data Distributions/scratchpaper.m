x = eyedat{cndlop}(1,:)*24+400;
y = eyedat{cndlop}(2,:)*24+300;
maxfixv = max(va(:,1))*2;
maxfixa = max(va(:,1))*2;
xpast = NaN(1,10);
ypast = NaN(1,10);
fixations = cell(1,length(eyedat));
time = 0;
extra = 0;
lag = 0;
for cndlop = 1:2:length(eyedat)
    type = NaN(1,length(x)); %1 fixation 0 sacccade
    for tt = 1:length(x);
        tic
        time = time + 0.005;
        xpast = [x(round(time/0.005)) xpast(1:end-1) ];
        ypast = [y(round(time/0.005)) ypast(1:end-1) ];
        if time > 0.05;
            velx = diff(xpast);
            vely = diff(ypast);
            vel = sqrt(velx.^2+vely.^2);
            accel = abs(diff(vel));
            accel = accel/.005^2/500;
            vel = vel/.005/5;
            velnow = round(mean(vel(end-7:end)));
            accelnow = round(mean(accel(end-7:end)));
            velnow(velnow > maxfixv) = maxfixv;
            accelnow(accelnow > maxfixa) = maxfixa;
            if full(statespace2D(velnow,accelnow)) == 1
                type(tt) = 1;
            else
                type(tt) = 0;
            end
        end
        tpass = toc;
        if tpass < 0.005
            extra = extra + 0.005-tpass;
        else
            lag = lag+0.005-tpass;
        end
    end
end

type(isnan(type)) = [];
saccind = find(type == 0);
sacgap = find(diff(saccind) > 1);
badsaccind = [];
for i = 1:length(sacgap)
    if i == 1;
        temp = saccind(1:sacgap(i));
        if length(temp) < 5
            badsaccind = [badsaccind 1:sacgap(i)];
        end
    elseif i == length(sacgap);
        temp = saccind(sacgap(i)+1:length(saccind));
        if length(temp) < 5
            badsaccind = [badsaccind sacgap(i)+1:length(saccind)];
        end
    else
        temp = saccind(sacgap(i-1)+1:sacgap(i));
        if length(temp) < 5
            badsaccind = [badsaccind sacgap(i-1)+1:sacgap(i)];
        end
    end
end
saccind(badsaccind) = [];
fixind = 10:length(x);
[~, ia, ~] = intersect(fixind,saccind);
fixind(ia) = [];
fixgap = find(diff(fixind) > 1);
badfixind = [];
for i = 1:length(fixgap)
    if i == 1;
        temp = fixind(1:fixgap(i));
        if length(temp) < 10
            badfixind = [badfixind 1:fixgap(i)];
        end
    elseif i == length(fixgap);
        temp = fixind(fixgap(i)+1:length(fixind));
        if length(temp) < 10
            badfixind = [badfixind fixgap(i)+1:length(fixind)];
        end
    else
        temp = fixind(fixgap(i-1)+1:fixgap(i));
        if length(temp) < 10
            badfixind = [badfixind fixgap(i-1)+1:fixgap(i)];
        end
    end
end
fixind(badfixind) = [];
length(find(diff(fixind) > 1))
toc
figure
plot(x,y,'g');
hold on
fixgap = find(diff(fixind) > 1);
for i = 1:length(fixgap)
    if i == 1;
        temp = fixind(1:fixgap(i));
        
    elseif i == length(fixgap);
        temp = fixind(fixgap(i)+1:length(fixind));
        
    else
        temp = fixind(fixgap(i-1)+1:fixgap(i));
    end
    plot(x(temp),y(temp),'r')
end
hold off