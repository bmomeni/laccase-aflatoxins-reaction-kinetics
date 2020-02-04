clear

f1 = 56000;
f2 = 54;
f3 = 32.8;
% FS = f3 + f1*TS./(TS+f2);
% Teq = f2*(F-f3)./(f1-F+f3);

%% 25U/ml
load('Laccase25U AFG2 and AFB1 kinetics experiment.mat')
c = 9; %column
r = 1:2; %row
SBG = mean(shiftdim(FL2(r,c,1:Nr),2),2);
% figure
% plot(SBG)

dt = 10/60; % (time interval, ten minutes)
E = 25; % enzyme conc. (U/mL)

% 30 ug/ml
c = 1:4; %column
r = 4; %row
S = shiftdim (FL2(r,c,1:Nr),1)-ones(4,1)*SBG';
S30 = mean(mean(S(:,1:2)));
S30t = S;
Tox30c = f2*(S-f3)./(f1-S+f3);
for cp = 1:50
    x = floor((Nr-31)*rand(1))+26;
    for cnt = 1:4
        T30(cp,cnt) = dt*x;
        Tox30(cp,cnt) = 30/S30*S(cnt,x);
        Tox30s(cp,cnt) = Tox30c(cnt,x);
        pf = polyfit(dt*(x:x+5),Tox30c(cnt,x:x+5),1);
        DegRate30(cp,cnt) = -1/E*pf(1);
    end
end

% 50 ug/ml
c = 1:4; %column
r = 5; %row1
S = shiftdim (FL2(r,c,1:Nr),1)-ones(4,1)*SBG';
S50 = mean(mean(S(:,1:2)));
S50t = S;
Tox50c = f2*(S-f3)./(f1-S+f3);
for cp = 1:100
    x = floor((Nr-6)*rand(1))+1;
    for cnt = 1:4
        T50(cp,cnt) = dt*x;
        Tox50(cp,cnt) = 50/S50*S(cnt,x);
        Tox50s(cp,cnt) = Tox50c(cnt,x);
        pf = polyfit(dt*(x:x+5),Tox50c(cnt,x:x+5),1);
        DegRate50(cp,cnt) = -1/E*pf(1);
    end
end

% 100 ug/ml
c = 1:4; %column
r = 6; %row1
S = shiftdim (FL2(r,c,1:Nr),1)-ones(4,1)*SBG';
S100 = mean(mean(S(:,1:2)));
S100t = S;
Tox100c = f2*(S-f3)./(f1-S+f3);
for cp = 1:100
    x = floor((Nr-6)*rand(1))+1;
    for cnt = 1:4
        T100(cp,cnt) = dt*x;
        Tox100(cp,cnt) = 100/S100*S(cnt,x);
        Tox100s(cp,cnt) = Tox100c(cnt,x);
        pf = polyfit(dt*(x:x+5),Tox100c(cnt,x:x+5),1);
        DegRate100(cp,cnt) = -1/E*pf(1);
    end
end

figure
plot(0:dt:96,Tox30c)
hold on
plot(0:dt:96,Tox50c)
plot(0:dt:96,Tox100c)
xlabel('Time (hrs)')
ylabel('AFG2 Conc. (\mug/ml)')

% KT = 500;
% DRM = 1.6;
eta = 2.7e-3;
TS = linspace(0,120,200);
% DegRateModel = DRM*TS./(TS+KT);
DegRateModel = eta*TS;
figure
plot(Tox30s,DegRate30,'.','color',[0.1 0.7 0.2])
hold on
plot(Tox50s,DegRate50,'.','color',[0.1 0.7 0.2])
% plot(Tox100s,DegRate100,'.','color',[0.1 0.7 0.2])
plot(TS,DegRateModel,'k:')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('Degradation rate (\mugU^-^1h^-^1)')
% title('Conc. from local normalization')
xlim([0 80])
ylim([0 0.25])

figure
plot(20*dt:dt:96,mean(Tox30c(:,21:Nr),1),'color',[0.1 0.7 0.2])
hold on
plot(20*dt:dt:96,mean(Tox50c(:,21:Nr),1),'color',[0.1 0.7 0.2])
plot(20*dt:dt:96,mean(Tox100c(:,21:Nr),1),'color',[0.1 0.7 0.2])
xlabel('Time (hrs)')
ylabel('AFB1 conc. (\mug/ml)')
ylim([0 80])

figure
plot(24*dt:dt:96,Tox30c(:,25:Nr),'color',[0.1 0.7 0.2])
hold on
% plot(20*dt:dt:96,mean(Tox50c(:,21:Nr),1),'color',[0.1 0.7 0.2])
plot(24*dt:dt:96,Tox100c(:,25:Nr),'color',[0.1 0.7 0.2])
xlabel('Time (hrs)')
ylabel('AFB1 conc. (\mug/ml)')
ylim([0 80])

figure
plot(T30,Tox30s,'.')
hold on
plot(T50,Tox50s,'.')
plot(T100,Tox100s,'.')
xlabel('Time (hrs)')
ylabel('AFG2 Conc. (\mug/ml)')

figure
plot(T30,DegRate30,'r.')
hold on
plot(T50,DegRate50,'g.')
plot(T100,DegRate100,'b.')
xlabel('Time (hrs)')
ylabel('Degradation rate (\mugU^-^1h^-^1)')
ylim([0 0.2])


