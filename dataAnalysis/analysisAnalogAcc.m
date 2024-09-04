function analysisAccelerometer(AVPlotDataObj)

%zeroGBias=[507384.049419002,489856.736624713,463392.702104808]';
%sensitivity=[34291.9851784505,35144.5526086612,35216.4584260847]';

zeroGBias=[472971.375013507,494947.777368294,493554.966062353]';
sensitivity=[34735.1293215227,34954.7241622606,34589.4195371905]';

%convert accelerometer measurements from voltage to acceleration: a=(V-zero_bias)/sensitivity
AVPlotDataObj.A(:) = bsxfun(@rdivide, bsxfun(@minus,AVPlotDataObj.A,zeroGBias),sensitivity);

%p = poseplot(eye(3),[5 5 5],ScaleFactor=0.5)

%figure;plot(squeeze(AVPlotDataObj.A)')
%AVPlotDataObj.A=AVPlotDataObj.A.*0.05481;

%figure;plot(sqrt(squeeze(AVPlotDataObj.A(1,1,:).^2+AVPlotDataObj.A(2,1,:).^2+AVPlotDataObj.A(3,1,:).^2)))
%fprintf('Total acceleration is %f.',mean(sqrt(squeeze(AVPlotDataObj.A(1,1,1:20000).^2+AVPlotDataObj.A(2,1,1:20000).^2+AVPlotDataObj.A(3,1,1:20000).^2))));
%set first sample to zero before start integrating
%AVPlotDataObj.A = bsxfun(@rdivide, bsxfun(@minus,AVPlotDataObj.A,mean(AVPlotDataObj.A(:,:,1:100),3)),sensitivity);

%AVPlotDataObj.A=cumtrapz(AVPlotDataObj.A,3); %Integrate acceleration to get velocity
%AVPlotDataObj.A=cumtrapz(AVPlotDataObj.A,3);

%{
"X": {
    "sensitivity": -3377.763012130533,
    "zero_g_bias_level": 490.0634399466667
  },
  "Y": {
    "sensitivity": -3386.7624159678026,
    "zero_g_bias_level": 463.3951879333333
  },
  "Z": {
    "sensitivity": 3297.7114988816197,
    "zero_g_bias_level": 507.3929007399999
  }


t = data(:,1);
acc = data(:,2);
acc = detrend(acc,%}'linear');
N = length(t);
dt = mean(diff(t));       % Average dt
fs = 1/dt;                  % Frequency [Hz] or sampling rate
% some additionnal high pass filtering
N = 4;
fc = 0.05; % Hz
[B,A] = butter(N,2*fc/fs,'high');
acc2 = filter(B,A,acc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocity = cumtrapz(dt,acc2);
velocity = detrend(velocity,'linear');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp     = cumtrapz(dt,velocity);




X_m=[           442899.07831234,
          469150.639233738,
          826769.064993549];

X_p=[480128.46884276,
          488016.032968798,
          828920.179462916]


Xsensitivity=(X_p-X_m)/2/9.80665;
Xbias=(X_p+X_m)/2;

aX=(AVPlotDataObj.A(1,1,:)-Xbias(1))/Xsensitivity(1);


folder='/media/sil2/Data/Nimrod/OE_headstage_calib_2024_02_07/headstage_calib_01_4918';
file={...
'X_negative_01_2024-02-07_10-32-17',  
'X_positive_01_2024-02-07_10-31-13', 
'Y_negative_01_2024-02-07_10-26-50', 
'Y_positive_01_2024-02-07_10-23-59',
'Z_negative_01_2024-02-07_10-20-06',  
'Z_positive_01_2024-02-07_10-18-35',...
};
ax=[1,1,2,2,3,3];
side=[1,2,1,2,1,2];
trans=[3,3,1,1,2,2];

folder='/media/sil2/Data/Nimrod/OE_headstage_calib_2024_02_07/headstage_calib_02_0817';
file={...
'X_negative_01_2024-02-07_10-37-15',
'X_positive_01_2024-02-07_10-36-44',
'Y_negative_01_2024-02-07_10-41-11',
'Y_positive_02_2024-02-07_10-39-30',
'Z_negative_01_2024-02-07_10-43-37',
'Z_positive_01_2024-02-07_10-43-02',
};

ax=[1,1,2,2,3,3];
side=[1,2,1,2,1,2];
trans=[1,1,2,2,3,3];

for i=1:numel(file)
    OE=OERecording([folder,filesep,file{i},filesep,'Record Node 102']);
    V=OE.getAnalogData([],0,OE.recordingDuration_ms);
    sensitivity(ax(i),side(i))=mean(V(trans(i),1,:));
end

sensi=-(sensitivity(:,1)-sensitivity(:,2))/2/9.80665;
bias=(sensitivity(:,1)+sensitivity(:,2))/2;
%}