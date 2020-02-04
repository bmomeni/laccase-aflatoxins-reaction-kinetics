clear

% TS = linspace(0,150,200);
f1 = 49800;
f2 = 300;
f3 = 268.0667;
% FS = f3 + f1*TS./(TS+f2);

%% 25U/ml
load('Laccase25U AFG2 and AFB1 kinetics experiment.mat')
c = 9; %column
r = 1:2; %row
SBG = mean(shiftdim(FL1(r,c,1:Nr),2),2); % background (laccase+buffer+methanol)
dt = 10/60; % (time interval, ten minutes)
E = 25; % enzyme conc. (U/mL)
pe = 50;

load('Laccase25U AFG2 and AFB1 kinetics experiment.mat')
% 30 ug/ml, laccase, toxin, buffer, and methanol
c = 2:4; %column
r = 1; %row1
% S = shiftdim (FL1(r,c,1:Nr),1)-ones(3,1)*SBG';
S = shiftdim (FL1(r,c,1:Nr),1)-SBG(1);
S30t = S;
S30 = mean(mean(S(:,1:2)));
Tox30c = f2*(S-f3)./(f1-S+f3);
for cp = 1:50
    x = floor((Nr-pe-1)*rand(1))+1;
    for cnt = 1:3
        T30(cp,cnt) = dt*x;
%         Tox30(cp,cnt) = 0.75*30/S30*S(cnt,x);
        Tox30s(cp,cnt) = Tox30c(cnt,x);
        pf = polyfit(dt*(x:x+pe),Tox30c(cnt,x:x+pe),1);
        DegRate30(cp,cnt) = -1/E*pf(1);
    end
end

% 50 ug/ml
c = 2:4; %column
r = 2; %row1
S = shiftdim (FL1(r,c,1:Nr),1)-ones(3,1)*SBG';
% S = shiftdim (FL1(r,c,1:Nr),1)-SBG(1);
S50t = S;
S50 = mean(mean(S(:,1:2)));
Tox50c = f2*(S-f3)./(f1-S+f3);
for cp = 1:50
    x = floor((Nr-pe-1)*rand(1))+1;
    for cnt = 1:3
        T50(cp,cnt) = dt*x;
%         Tox50(cp,cnt) = 0.75*50/S50*S(cnt,x);
        Tox50s(cp,cnt) = Tox50c(cnt,x);
        pf = polyfit(dt*(x:x+pe),Tox50c(cnt,x:x+pe),1);
        DegRate50(cp,cnt) = -1/E*pf(1);
    end
end

% 100 ug/ml
c = [2 4 5]; %column
r = 3; %row1
S = shiftdim (FL1(r,c,1:Nr),1)-ones(3,1)*SBG';
% S = shiftdim (FL1(r,c,1:Nr),1)-SBG(1);
S100t = S;
S100 = mean(mean(S(:,1:2)));
Tox100c = f2*(S-f3)./(f1-S+f3);
for cp = 1:50
    x = floor((Nr-pe-1)*rand(1))+1;
    for cnt = 1:3
        pf = polyfit(dt*(x:x+pe),S(cnt,x:x+pe),1);
        T100(cp,cnt) = dt*x;
%         Tox100(cp,cnt) = 0.64*100/S100*S(cnt,x);
        Tox100s(cp,cnt) = Tox100c(cnt,x);
        pf = polyfit(dt*(x:x+pe),Tox100c(cnt,x:x+pe),1);
        DegRate100(cp,cnt) = -1/E*pf(1);
    end
end

figure
plot(25*dt:dt:96,mean(Tox30c(:,26:Nr),1),'color',[0 0.2 1])
hold on
plot(25*dt:dt:96,mean(Tox50c(:,26:Nr),1),'color',[0 0.2 1])
plot(25*dt:dt:96,mean(Tox100c(:,26:Nr),1),'color',[0 0.2 1])
xlabel('Time (hrs)')
ylabel('AFB1 conc. (\mug/ml)')
ylim([0 140])

% figure
% plot(0:dt:96,30/S30*S30t,'color',[0 0.2 1])
% hold on
% plot(0:dt:96,50/S50*S50t,'color',[0 0.2 1])
% plot(0:dt:96,100/S100*S100t,'color',[0 0.2 1])
% xlabel('Time (hrs)')
% ylabel('AFB1 conc. (\mug/ml)')
% ylim([0 120])

% figure
% for ccnt = 1:3
% plot(0:dt:96,0.75*interp1(FS,TS,S30t(ccnt,:)),'color',[0 0.2 1])
% hold on
% plot(0:dt:96,0.75*interp1(FS,TS,S50t(ccnt,:)),'color',[0 0.2 1])
% plot(0:dt:96,0.64*interp1(FS,TS,S100t(ccnt,:)),'color',[0 0.2 1])
% end
% xlabel('Time (hrs)')
% ylabel('AFB1 conc. (\mug/ml)')
% ylim([0 100])

% KT = 500;
% DRM = 0.18;
eta = 3.0e-4;
TS = linspace(0,150,200);
% DegRateModel = DRM*TS./(TS+KT);
DegRateModel = eta*TS;
figure
plot(Tox30s,DegRate30,'.','color',[0 0.2 1])
hold on
plot(Tox50s,DegRate50,'.','color',[0 0.2 1])
plot(Tox100s,DegRate100,'.','color',[0 0.2 1])
plot(TS,DegRateModel,'k:')
xlabel('AFB1 Conc. (\mug/ml)')
ylabel('Degradation rate (\mugU^-^1h^-^1)')
xlim([0 150])
ylim([0 0.05])

figure
plot(T30,DegRate30,'.','color',[0.4 0.4 0.4])
hold on
plot(T50,DegRate50,'.','color',[0.4 0.4 0.4])
plot(T100,DegRate100,'.','color',[0.4 0.4 0.4])
xlabel('Time (hrs)')
ylabel('Degradation rate (\mugU^-^1h^-^1)')


% %% 10U/ml
% load('Laccase10U AFG2 and AFB1 kinetics experiment.mat')
% c = 9; %column
% r = 1:6; %row
% SBG = mean(shiftdim(FL1(r,c,1:Nr),2),2);
% % figure
% % plot(SBG)
% 
% load('Laccase_AFB1_kinetics_experiment.mat')
% dt = 10/60; % (time interval, ten minutes)
% E = 10; % enzyme conc. (U/mL)
% 
% % 3 ug/ml
% c = 4:6; %column
% r = 3; %row
% S1 = shiftdim (FL1(r,c,1:Nr),1)-ones(3,1)*SBG';
% 
% for cnt = 1:3
%     pf = polyfit(dt*(1:50),S1(cnt,25:74),1);
%     DegRate(1,cnt) = -1/E*pf(1);
% end
% 
% % 30 ug/ml
% c = 1:3; %column
% r = 2; %row
% S2 = shiftdim (FL1(r,c,1:Nr),1)-ones(3,1)*SBG';
% 
% for cnt = 1:3
%     pf = polyfit(dt*(1:50),S2(cnt,25:74),1);
%     DegRate(2,cnt) = -1/E*pf(1);
% end
% 
% load('Laccase10U AFG2 and AFB1 kinetics experiment.mat')
% % 30 ug/ml
% c = 1:5; %column
% r = 1; %row1
% S5 = shiftdim (FL1(r,c,1:Nr),1)-ones(5,1)*SBG';
% 
% for cnt = 1:5
%     pf = polyfit(dt*(1:50),S5(cnt,25:74),1);
%     DegRate(5,cnt) = -1/E*pf(1);
% end
% 
% % 50 ug/ml
% c = 1:5; %column
% r = 2; %row1
% S3 = shiftdim (FL1(r,c,1:Nr),1)-ones(5,1)*SBG';
% 
% for cnt = 1:5
%     pf = polyfit(dt*(1:50),S3(cnt,25:74),1);
%     DegRate(3,cnt) = -1/E*pf(1);
% end
% 
% % 100 ug/ml
% c = 1:5; %column
% r = 3; %row1
% S4 = shiftdim (FL1(r,c,1:Nr),1)-ones(5,1)*SBG';
% 
% for cnt = 1:5
%     pf = polyfit(dt*(1:50),S4(cnt,25:74),1);
%     DegRate(4,cnt) = -1/E*pf(1);
% end
% 
% figure
% plot(0:dt:96,S1)
% hold on
% plot(0:dt:96,S2)
% plot(0:dt:96,S3)
% plot(0:dt:96,S4)
% 
% KT = 22;
% DRM = 3.6;
% TS = linspace(0,120,200);
% DegRateModel = DRM*TS./(TS+KT);
% T1 = [30 50 100];
% T2 = [3  30];
% figure
% errorbar(T1,mean(DegRate(2:4,:)'),std(DegRate(2:4,:)'),'bo')
% hold on
% errorbar(T2,mean(DegRate([1 5],1:3)'),std(DegRate([1 5],1:3)'),'bo')
% plot(TS,DegRateModel,'k:')
% xlabel('AFB1 Conc. (\mug/ml)')
% ylabel('Degradation rate (a.u.)')
