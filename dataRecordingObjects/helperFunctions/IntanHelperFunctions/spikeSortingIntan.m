function [fileNames]=spikeSortingMultiChannel(dataRecordingObj,varargin)
%default variables
overwriteSpikeExtraction=0;
overwriteFeatureExtraction=0;
overwriteClustering=0;
overwriteMerging=1;
overwriteFitting=1;
overwritePostProcessignAnalysis=0;

upSamplingFrequency=60000;
selectedChannels=[];

gridSize=5; %must be an odd number. If empty extracts using all channels
gridExt=2;

maxSpikesToCluster=10000; %maximum number of spikes to feed to feature extraction and clustering (not to template matching)

%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
    end
    return;
end

[~, name, ~] = fileparts(dataRecordingObj.recordingName{1});
sortingDir=[dataRecordingObj.recordingDir filesep name '_spikeSort'];

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

%make directory in recording folder with spike sorting data
if ~exist(sortingDir,'dir')
    mkdir(sortingDir);
end

%determine channels
if isempty(selectedChannels)
    selectedChannels=dataRecordingObj.channelNumbers; %take all channel in the recording
end
nCh=numel(selectedChannels);
chPar.s2r=selectedChannels; %transformation between the serial channel number and the real ch number
chPar.r2s(selectedChannels)=1:nCh; %transformation between the real channel number and the serial ch number

%get channel layout
chPar.En=dataRecordingObj.chLayoutNumbers;

%limit layout only to the serial channels selected by user
chPar.rEn=chPar.En;
chPar.rEn(~isnan(chPar.En))=chPar.r2s(chPar.En(~isnan(chPar.En)));
chPar.rEn(chPar.rEn==0)=NaN;
[nRowsTmp,nColsTmp]=size(chPar.En);

gridSizeExt=gridSize+gridExt*2;

%initiate arrays
chPar.surChExt=cell(1,nCh);chPar.pValidSurChExt=cell(1,nCh);chPar.surChExtVec=cell(1,nCh);chPar.pCenterCh=zeros(1,nCh);chPar.nValidChExt=zeros(1,nCh);chPar.pSurCh=cell(1,nCh);chPar.pSurChOverlap=cell(1,nCh);
if ~isempty(gridSizeExt)
    gridExtension=(gridSize-1)/2;
    overheadGridExtension=(gridSizeExt-1)/2;
    EnExt=NaN(nRowsTmp+overheadGridExtension*2,nColsTmp+overheadGridExtension*2);
    EnExt(1+overheadGridExtension:end-overheadGridExtension,1+overheadGridExtension:end-overheadGridExtension)=chPar.rEn;
    for i=1:nCh
        [x,y]=find(EnExt==i);
        
        %find the surrounding channels on which feature extraction will be performed
        surCh=EnExt(x-gridExtension:x+gridExtension,y);  %this is different for intan because I don't want to sort across electrodes
        pValidSurCh=find(~isnan(surCh)); %do not remove find
        
        %find the channels that are overhead step from the central channel - these are the channels who's waveforms should be checked for merging
        surChOverlap=EnExt(x-gridExt:x+gridExt,y);
        surChOverlap(gridExt+1)=NaN;  %is this removal of central channel?
        pValidSurChOverlap=find(~isnan(surChOverlap)); %do not remove find
        
        %find the extended channels for merging of the same neurons detected on nearby channels
        chPar.surChExt{i}=EnExt(x-overheadGridExtension:x+overheadGridExtension,y);
        chPar.pValidSurChExt{i}=find(~isnan(chPar.surChExt{i})); %do not remove find
        chPar.nValidChExt(i)=numel(chPar.pValidSurChExt{i});
        
        chPar.surChExtVec{i}=chPar.surChExt{i}(chPar.pValidSurChExt{i}(:))';
        chPar.pCenterCh(i)=find(chPar.surChExtVec{i}==i); %the position of the central channel in surChExtVec
        [~,chPar.pSurCh{i}]=intersect(chPar.surChExtVec{i},surCh(pValidSurCh(:)));
        [~,chPar.pSurChOverlap{i}]=intersect(chPar.surChExtVec{i},surChOverlap(pValidSurChOverlap(:)));        
    end
end
%map channel intersections
for i=1:nCh
    for j=chPar.surChExtVec{i}(chPar.pSurChOverlap{i}) %the trivial case : (i,i) is also included (can be remove in not required)
        [chPar.sharedChNames{i}{j},chPar.pSharedCh1{i}{j},chPar.pSharedCh2{i}{j}]=intersect(chPar.surChExtVec{i},chPar.surChExtVec{j});
    end
end

%initiate variables
t=[];ic=[];avgWaveform=[];

%check conditions for recalculating the different stages of spike sorting. Notice that for the 3 first procedures only, calcululation of a subset of uncalculated channels is possiblee
for i=1:nCh
    fileNames.spikeDetectionFile{i}=[sortingDir filesep 'ch_' num2str(chPar.s2r(i)) '_spikeDetection.mat'];
    fileNames.spikeDetectionExist(i)=exist(fileNames.spikeDetectionFile{i},'file');
    
    fileNames.featureExtractionFile{i}=[sortingDir filesep 'ch_' num2str(chPar.s2r(i)) '_featureExtraction.mat'];
    fileNames.featureExtractionExist(i)=exist(fileNames.featureExtractionFile{i},'file');
    
    fileNames.clusteringFile{i}=[sortingDir filesep 'ch_' num2str(chPar.s2r(i)) '_clustering.mat'];
    fileNames.clusteringExist(i)=exist(fileNames.clusteringFile{i},'file');
end

fileNames.avgWaveformFile=[sortingDir filesep 'AllClusteredWaveforms.mat'];

fileNames.mergedAvgWaveformFile=[sortingDir filesep 'AllMergedWaveforms.mat'];
fileNames.mergedAvgWaveformExist=exist(fileNames.mergedAvgWaveformFile,'file');

fileNames.fittingFile=[sortingDir filesep 'spikeSorting.mat'];
fileNames.fittingExist=exist(fileNames.fittingFile,'file');

fileNames.postProcessignAnalysisFile=[sortingDir filesep 'postProcessignAnalysis.mat'];
fileNames.postProcessignAnalysisExist=exist(fileNames.postProcessignAnalysisFile,'file');

%the overwriting of files needed, treat all files as non existing
if overwriteSpikeExtraction,fileNames.spikeDetectionExist(:)=0;end;
if overwriteFeatureExtraction,fileNames.featureExtractionExist(:)=0;end;
if overwriteClustering,fileNames.clusteringExist(:)=0;end;
if overwriteMerging,fileNames.mergedAvgWaveformExist=0;end;
if overwriteFitting,fileNames.fittingExist=0;end;
if overwritePostProcessignAnalysis,fileNames.postProcessignAnalysisExist=0;end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  spike Detection  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
if all(fileNames.spikeDetectionExist)  %check for the existence of spike shapes
    disp('Sorting will be preformed on previously detected waveforms');
else
    spikeDetectionNSK(dataRecordingObj,chPar,fileNames,'upSamplingFrequency',upSamplingFrequency,'gridSize',gridSize,'gridExt',gridExt,...
        'dataFileFolder',sortingDir,varargin{:});
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  feature extraction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
if all(fileNames.featureExtractionExist)  %check for the existence of spike shapes
    disp('Sorting will be preformed on previously extracted features');
else
    spikeFeatureExtractionNSK(chPar,fileNames,'sortingDir',sortingDir,'maxSpikesToCluster',maxSpikesToCluster,varargin{:});
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  clustering  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if all(fileNames.clusteringExist)  %check for the existence of spike shapes
    disp('Sorting will be preformed on previously extracted clusters');
else
    spikeClusteringNSK(dataRecordingObj,chPar,fileNames,'sortingDir',sortingDir,'upSamplingFrequency',upSamplingFrequency,'maxSpikesToCluster',maxSpikesToCluster,varargin{:});
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Merging duplicate neurons  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if fileNames.mergedAvgWaveformExist  %check for the existence of spike shapes
    disp('Sorting will be preformed on previously merged clusters');
else
    spikeMergingClustersNSK(chPar,fileNames,'gridSize',gridSize,'gridExt',gridExt,'sortingDir',sortingDir);
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Fitting duplicate neurons  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if fileNames.fittingExist  %check for the existence of spike shapes
    disp('No fitting performed!!!');
else
    spikeFittingNSK(chPar,fileNames,'gridSize',gridSize,'gridExt',gridExt,'sortingDir',sortingDir,'upSamplingFrequency',upSamplingFrequency);
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  General final plots and access sorting quality  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if fileNames.postProcessignAnalysisExist  %check for the existence of spike shapes
    disp('No post analysis performed!!!');
else
    spikePostProcessingNSK(dataRecordingObj,chPar,fileNames,'sortingDir',sortingDir,'upSamplingFrequency',upSamplingFrequency);
end
toc;
