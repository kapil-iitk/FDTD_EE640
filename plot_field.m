clc;
clear all;
EData=csvread('EData.csv');
xx=1:201;;

for ii=1:1000
yy=EData(ii,:);
clf;
plot(xx,yy);
axis([0 251 -1 1])
pause(0.01);
end