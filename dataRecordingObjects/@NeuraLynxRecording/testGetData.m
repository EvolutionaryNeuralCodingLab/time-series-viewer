recName='D:\PostDoc\Experiments\test example\CSC2.ncs';
[timeStampsData, dataSamples] = Nlx2MatCSC(recName,[1 0 0 0 1],0,1,[]);
dataSamples=dataSamples(:);
t_ms=(1:numel(dataSamples))/32000*1000;

test=NeuraLynxRecording('D:\PostDoc\Experiments\test example');
dataSamples=dataSamples*test.A2DBitVolts*1e6;
win_ms=32;
winSamples=win_ms/1000*32000;

figure;
%Recording start
[V,t]=test.getData(2,[-16 0],win_ms);
subplot(2,2,1);
plot(t_ms(1:winSamples),dataSamples(1:winSamples));hold on;
plot(t-16,squeeze(V(1,1,:))','or');
plot(t-0,squeeze(V(1,2,:))','sg');

%Recording middle
[V2,t2]=test.getData(2,[1000 1010],win_ms);
p=find(t_ms>1000,1,'first');
subplot(2,2,2);
plot(t_ms(p:(p+winSamples-1)),dataSamples(p:(p+winSamples-1)));hold on;
plot(t2+1000,squeeze(V2(1,1,:))','or');
plot(t2+1010,squeeze(V2(1,2,:))','sg');

%Recording end
tStart_ms=t_ms(end-winSamples+1);
[V3,t3]=test.getData(2,[tStart_ms-32 tStart_ms tStart_ms+16],win_ms);
subplot(2,2,3);
plot(t_ms((end-winSamples+1):end),dataSamples((end-winSamples+1):end));hold on;
plot(t3+tStart_ms-32,squeeze(V3(1,1,:))','*k');
plot(t3+tStart_ms,squeeze(V3(1,2,:))','or');
plot(t3+tStart_ms+16,squeeze(V3(1,3,:))','sg');