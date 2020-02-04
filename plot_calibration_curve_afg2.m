clear

T = [100 50 25 12.5 6.25 3.125 1.6125 0.8625 0.43125 0];
Ns = 5;
%% 25U/ml
load('AFG2 and AFB1 concentration curve.mat')
c = 1:3; %column
r = 5; %row
G2(:,1) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
c = 7:9; %column
r = 5; %row
G2(:,2) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
G2(1,3) = mean(shiftdim(FL2(5,10,1:Ns),1),2);
c = 1:2; %column
r = 6; %row
G2(2:3,3) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
c = 3:5; %column
r = 6; %row
G2(:,4) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
c = 6:8; %column
r = 6; %row
G2(:,5) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
c = 9:10; %column
r = 6; %row
G2(1:2,6) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
G2(3,6) = mean(shiftdim(FL2(7,1,1:Ns),1),2);
c = 2:4; %column
r = 7; %row
G2(:,7) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
c = 5:7; %column
r = 7; %row
G2(:,8) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
c = 8:10; %column
r = 7; %row
G2(:,9) = mean(shiftdim(FL2(r,c,1:Ns),1),2);
c = 1:3; %column
r = 1; %row
G2(:,10) = mean(shiftdim(FL2(r,c,1:Ns),1),2);

G2(3,2) = NaN;
G2(3,3) = NaN;

TS = linspace(0,110,200);
f1 = 60000;
f2 = 54;
f3 = mean(G2(:,10));
FS = f3 + f1*TS./(TS+f2);

figure
plot(TS,0.001*FS,'r:')
hold on
errorbar(T,0.001*nanmean(G2),0.001*nanstd(G2),'bo')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('RFU (gain 50)')

figure
plot(T,0.001*G2,'bo')
xlabel('AFG2 Conc. (\mug/ml)')
ylabel('FL2 reading')

figure
plot(1/6*(1:80),shiftdim(FL2(2,1:3,1:80),1))
xlabel('Time (hrs)')
ylabel('FL2 (gain 50)')


