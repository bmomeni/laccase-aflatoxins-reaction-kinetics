clear

T = [100 50 25 12.5 6.25 3.125 1.6125 0.8625 0.43125 0];
Ns = 5;
%% 25U/ml
load('AFG2 and AFB1 concentration curve.mat')
c = 1:3; %column
r = 2; %row
B1(:,1) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
c = 7:9; %column
r = 2; %row
B1(:,2) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
B1(1,3) = mean(shiftdim(FL1(2,10,1:Ns),1),2);
c = 1:2; %column
r = 3; %row
B1(2:3,3) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
c = 3:5; %column
r = 3; %row
B1(:,4) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
c = 6:8; %column
r = 3; %row
B1(:,5) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
c = 9:10; %column
r = 3; %row
B1(1:2,6) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
B1(3,6) = mean(shiftdim(FL1(4,1,1:Ns),1),2);
c = 2:4; %column
r = 4; %row
B1(:,7) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
c = 5:7; %column
r = 4; %row
B1(:,8) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
c = 8:10; %column
r = 4; %row
B1(:,9) = mean(shiftdim(FL1(r,c,1:Ns),1),2);
c = 1:3; %column
r = 1; %row
B1(:,10) = mean(shiftdim(FL1(r,c,1:Ns),1),2);

TS = linspace(0,110,200);
f1 = 51700;
f2 = 300;
f3 = mean(B1(:,10));
FS = f3 + f1*TS./(TS+f2);

figure
plot(TS,0.001*FS,'r:')
hold on
errorbar(T,0.001*nanmean(B1),0.001*nanstd(B1),'bo')
xlabel('AFB1 Conc. (\mug/ml)')
ylabel('RFU (gain 65)')

figure
plot(T,0.001*B1,'bo')
xlabel('AFB1 Conc. (\mug/ml)')
ylabel('RFU (gain 65)')


