classdef filterData
    properties
        samplingFrequency=[]; %Sampling rate
        filterDimension=1; %the dimensionality of the filter
        filterOrder=8; %filter order
        filterDesign='butter'; %filter type 'ellip'/'butter'/'kaiserwin'
        rippleInPassband=0.5; %Ap
        attenuationInHighpass=4; %Ast
        attenuationInLowpass=4; %Ast
        downSamplingFactor=[]; %down sampling factor for decimation
        lowPassStopCutoff=[]; %low-pass cutoff frequency
        lowPassPassCutoff=[];
        highPassStopCutoff=[]; %high-pass cutoff frequency
        highPassPassCutoff=[];
        highPassCutoff=[];%low-pass cutoff frequency - for use with butterworth
        lowPassCutoff=[];%high-pass cutoff frequency - for use with butterworth
        upSamplingFactor=[]; %up sampling factor for decimation
        currentDesignedFilter=[];
        filteredSamplingFrequency=[];
        padding=true;
        padWithMirrorTraces=true;
        SOS_b=[];
        G_a=[];
        vars=[];
    end
    properties (Hidden)
        dF=[] %fdesign filter specification object
        F=[] %design filter object
    end
    methods
        function nPaddingSamples=getNumberOfRequiredPaddingSamples(obj)
            nPaddingSamples=impzlength(obj.SOS_b);
        end
        function [obj]=designHighPass(obj)
            %design highpass filter
            %Usage: [HP]=HP.designHighPass(cutoff)
            %Input : cutoff - the cutoff frequency [Hz]
            %Output: HP - the highpass filter object
            obj.currentDesignedFilter='highpass';
            if ~isempty(obj.highPassPassCutoff)
                obj.vars={'Fst,Fp,Ast,Ap',obj.highPassStopCutoff,obj.highPassPassCutoff,obj.attenuationInHighpass,obj.rippleInPassband};
            elseif ~isempty(obj.highPassCutoff)
                obj.vars={'N,F3db',obj.filterOrder,obj.highPassCutoff};
                obj.filterDesign='butter';
            end
            obj=obj.runDesign;
        end
        function [obj]=designLowPass(obj)
            %design lowpass filter
            %Usage: [LP]=LP.designLowPass(cutoff)
            %Input : cutoff - the cutoff frequency [Hz]
            %Output: LP - the lowpass filter object
            if obj.lowPassCutoff>obj.samplingFrequency
                disp(['Error: low pass cutoff is above the sampling frequency, changing cutoff to Nyquist: ' num2str(obj.samplingFrequency/2)]);
                obj.lowPassCutoff=obj.samplingFrequency/2;
            end
            obj.currentDesignedFilter='lowpass';
            if ~isempty(obj.lowPassPassCutoff)
                obj.vars={'Fp,Fst,Ap,Ast',obj.lowPassPassCutoff,obj.lowPassStopCutoff,obj.rippleInPassband,obj.attenuationInLowpass};
            elseif ~isempty(obj.lowPassCutoff)
                obj.vars={'N,F3db',obj.filterOrder,obj.lowPassCutoff};
                obj.filterDesign='butter';
            end
            obj=obj.runDesign;
        end
        function [obj]=designBandPass(obj)
            %design bandpass filter
            %Usage: [BP]=BP.designBandPass(lowPassCutoff,highPassCutoff)
            %Input : highPassCutoff - the highpass cutoff frequency [Hz]
            %        lowPassCutoff - the lowpass cutoff frequency [Hz]
            %Output: BP - the bandpass filter object
            if ~isempty(obj.lowPassCutoff)
                if obj.lowPassCutoff>obj.samplingFrequency
                    disp(['Error: low pass cutoff is above the sampling frequency, changing cutoff to Nyquist: ' num2str(obj.samplingFrequency/2)]);
                    obj.lowPassCutoff=obj.samplingFrequency/2;
                elseif obj.highPassCutoff>=obj.lowPassCutoff
                    disp('Error: low pass cutoff is lower than highpass cutoff, please change settings');
                    return;
                end
            end
            obj.currentDesignedFilter='bandpass';
            if ~isempty(obj.lowPassPassCutoff)
                obj.vars={'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',obj.highPassStopCutoff,obj.highPassPassCutoff,obj.lowPassPassCutoff,obj.lowPassStopCutoff,obj.attenuationInHighpass,obj.rippleInPassband,obj.attenuationInLowpass};
            elseif ~isempty(obj.highPassCutoff)
                obj.vars={'N,F3dB1,F3dB2',obj.filterOrder,obj.highPassCutoff,obj.lowPassCutoff};
                obj.filterDesign='butter';
            end
            obj=obj.runDesign;
        end
        function [obj]=designBandStop(obj)
            %design bandstop filter
            %Usage: [BS]=BS.designBandStop(lowPassCutoff,highPassCutoff)
            %Input : highPassCutoff - the high cutoff frequency [Hz]
            %        lowPassCutoff - the low cutoff frequency [Hz]
            %Output: BS - the bandstop filter object
            obj.currentDesignedFilter='bandstop';
            if ~isempty(obj.lowPassPassCutoff)
                obj.vars={'Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2',obj.highPassPassCutoff,obj.highPassStopCutoff,obj.lowPassStopCutoff,obj.lowPassPassCutoff,obj.attenuationInHighpass,obj.rippleInPassband,obj.attenuationInLowpass};
            elseif ~isempty(obj.highPassCutoff)
                obj.vars={'N,F3dB1,F3dB2',obj.filterOrder,obj.highPassCutoff,obj.lowPassCutoff};
                obj.filterDesign='butter';
            end
            obj=obj.runDesign;
        end
        
        function [obj]=designNotch(obj,centerFreq,qualityFactor,assbandAttenuation)
            if ~strcmp(obj.filterDesign,'cheby1')
               disp(['Error: Design method should be cheby1 for the notch filter, add obj.filterDesign=''cheby1'';']);
               return;
            end
            obj.currentDesignedFilter='notch';
            if nargin==1
                qualityFactor=10;
                assbandAttenuation=1;
                centerFreq=50;
            else
                if nargin~=4, error('Input parameters missing: centerFreq,qualityFactor,assbandAttenuation');end
            end
            if round(obj.filterOrder/2)*2~=obj.filterOrder
                error('Notch must be of even order');
            end
            obj.filterOrder=2;
            disp('Notice that filter order is changed to 2 for notch to eliminate ripples');
            obj.vars={'N,F0,Q,Ap',obj.filterOrder,centerFreq,qualityFactor,assbandAttenuation};
            obj=obj.runDesign;
        end
        
        function [obj]=designUpSample(obj)
            if isempty(obj.upSamplingFactor)
                disp('up SamplingFactor property missing!');
                return;
            end

            obj.dF=fdesign.interpolator(obj.upSamplingFactor);
            %obj.dF = fdesign.polysrc(obj.upSamplingFactor,1,'Fractional Delay','Np',2,obj.samplingFrequency);
            
            obj.F=design(obj.dF);
            
            [obj.SOS_b,obj.G_a]=tf2sos(obj.F.numerator,1);
            
            obj.currentDesignedFilter='upSample';
            obj.filteredSamplingFrequency=obj.upSamplingFactor*obj.samplingFrequency;
            %y = interp(x,r);
        end
        function [obj]=designDownSample(obj)
            %design downsampling filter
            %Usage: [DS]=DS.designDownSample(downSamplingFactor)
            %Input : downSamplingFactor - the down sampling factor
            %Output: DS - the down sampling filter object
            %method should be rewritten with the filter design toolbox
            if isempty(obj.downSamplingFactor)
                disp('downSamplingFactor property missing!');
                return;
            end
            if isempty(obj.lowPassCutoff)
                obj.lowPassCutoff=0.8*(obj.samplingFrequency/2)/obj.downSamplingFactor;
                disp(['cutting frequency for downsampling was set to: ' num2str(obj.lowPassCutoff)]);
            end
            obj.filterOrder=8;
            [obj]=designLowPass(obj);
            obj.currentDesignedFilter='downSample';
            obj.filteredSamplingFrequency=obj.samplingFrequency/obj.downSamplingFactor;
            
            %obj.dF=fdesign.decimator(downSamplingFactor,'lowpass','N,F3db',obj.filterOrder,CutOffFreq,obj.samplingFrequency);
        end
        function filterDesigns=getDesigns(obj)
            filterDesigns = designmethods(obj.dF);
        end
        function [obj]=runDesign(obj)
            obj.dF=fdesign.(obj.currentDesignedFilter)(obj.vars{:},obj.samplingFrequency);
            obj.F=design(obj.dF,obj.filterDesign);
            obj.filterOrder=order(obj.F);
            obj.filteredSamplingFrequency=obj.samplingFrequency;
            obj.SOS_b=obj.F.sosMatrix;
            obj.G_a=obj.F.scaleValues;
        end
        function [filterMethods,filterNames]=getFilters(obj)
            allMethods=methods(class(obj));
            filterMethods=[allMethods(strncmp('design',allMethods,6))];
            filterNames=lower(cellfun(@(x) x(7:end),filterMethods,'UniformOutput',0)); %remove the "design" from the string
            filterMethods=[filterMethods;'GUIdesign'];
            filterNames=[filterNames;{'GUIdesign'}];
        end
        function [fData,T_ms]=getFilteredData(obj,Data)
            %filter data using the fil object
            %Usage: [fData,T_ms]=obj.getFilteredData(Data)
            %Input : Data - data to be filtered (format: [nChannels x nTrials x nSamples])
            %Output: fData - the filtered data
            [nCh,nTrial,nSamples]=size(Data);
            if obj.padding
                paddSamp=obj.getNumberOfRequiredPaddingSamples;
                if paddSamp>nSamples
                    disp('Trace is too short for padding, will padd with available samples, try also padding with zeros');
                    paddSamp=nSamples;
                end
                if obj.padWithMirrorTraces
                    Data=cat(3,Data(:,:,paddSamp:-1:1),Data,Data(:,:,end:-1:(end-paddSamp+1)));
                else
                    Data=cat(3,zeros(nCh,nTrial,paddSamp),Data,zeros(nCh,nTrial,paddSamp));
                end
                nSamples=nSamples+2*paddSamp;
            end
            if ~isempty(obj.F)
                fData=zeros(nCh,nTrial,nSamples);
                if nTrial>1 && nCh>1
                    Data=shiftdim(Data,2);
                    if nCh<=nTrial
                        for i=1:nCh
                            fData(i,:,:)=filtfilt(obj.SOS_b, obj.G_a, squeeze(Data(:,i,:)))';
                        end
                    else
                        for i=1:nTrial
                            fData(:,i,:)=filtfilt(obj.SOS_b, obj.G_a, squeeze(Data(:,:,i)))';
                        end
                    end
                elseif nCh==1
                    [fData(:)] = filtfilt(obj.SOS_b, obj.G_a, shiftdim(Data,1)')';
                else
                    [fData(:)] = filtfilt(obj.SOS_b, obj.G_a, shiftdim(Data,2))';
                end
                
                if obj.padding
                    nSamples=nSamples-2*paddSamp;
                    fData=fData(:,:,(paddSamp+1):(paddSamp+nSamples));
                end
                
                if ~isempty(obj.downSamplingFactor)
                    fData=fData(:,:,1:obj.downSamplingFactor:nSamples);
                end
                if nargout==2
                    T_ms=(1:size(fData,3))*(1e3/obj.filteredSamplingFrequency);
                end
            else
                disp('No filter was designed!!!!');
            end
        end
        function generateFilteredDataFromDataRecordingObject(obj,dataRecordingObj,fileName,startTime,window,channels)
            if nargin==5
                channels=[];
            elseif nargin==3
                channels=[];window=[];startTime=[];
            elseif nargin==2
                channels=[];window=[];startTime=[];fileName=[];
            elseif nargin==4 || nargin==1 || nargin>6
                error('The number of inputs was not valid!!!!');
            end
            if exist(fileName,'file') || exist([fileName '.mat'],'file')
                choice = questdlg('File name already exists, do you want to overwrite?');
                if length(choice)==3
                    delete(fileName);
                else
                    fileName=[];
                end
            end
            %open matlab file for writing and get pointer to file
            pointerToMatFile=matfile(fileName,'Writable',true);
            %copy all properties of dataRecording to file
            dataRecordingMetadata=?dataRecording;
            for i=1:numel(dataRecordingMetadata.PropertyList)
                tmpPropName=dataRecordingMetadata.PropertyList(i).Name;
                pointerToMatFile.(tmpPropName)=dataRecordingObj.(tmpPropName);
            end
            %add filter properties
            filterProps=properties(obj);
            for i=1:numel(filterProps)
                pointerToMatFile.(filterProps{i})=obj.(filterProps{i});
            end
            %examine input parameters
            if isempty(fileName)
                fileName = inputdlg({'full file name (including directory)'});
                if exist(fileName,'file')
                    error('The chosen fileName already exists');
                end
            end
            if isempty(channels)
                channels=dataRecordingObj.channelNumbers;
            end
            if isempty(startTime)
                startTime=0;
            end
            if isempty(window)
                window=dataRecordingObj.recordingDuration_ms;
            end
            %start filtering and saving data
            h = waitbar(0.5,'Saving triggers to file...');
            [t_ms]=dataRecordingObj.getTrigger();
            nTrig=numel(t_ms);
            for i=1:numel(t_ms)
                waitbar(i / nTrig)
                pointerToMatFile.(['Trigger' num2str(i)])=t_ms{i};
            end
            close(h);
            
            h = waitbar(0,'Saving filtered data to file...');
            nCh=numel(channels);
            for i=1:nCh
                waitbar(i / nCh)
                [V_uV,t_ms]=dataRecordingObj.getData(channels(i),startTime,startTime+window);
                pointerToMatFile.Data(i,1:numel(V_uV))=squeeze(obj.getFilteredData(V_uV))';
            end
            close(h);
        end
        function plotFilter(obj)
            %plot the designed filter
            %Usage: fil.plotFilter;
            %Input :
            %Output:
            if ~isempty(obj.F)
                fvtool(obj.F);
            else
                disp('No filter was designed');
            end
        end
        function [obj]=GUIdesign(obj)
            %use Matlab filter builder to build the object (the filter should be named F)
            %Usage: [fil]=fil.getFilteredData;
            %Input :
            %Output: fil - the filter obj
            h=msgbox({'Please remember to do the following for applying the filter:',...
                '1) Export to a filter object named F','2) Update the data sampling frequency'},'Attention','warn','modal');
            uiwait(h);
            
            myfilterbuilder;
            %h=fdatool;
            obj.F=evalin('base','F');
            obj.dF=obj.F.getfdesign;
            if isprop(obj.F,'Numerator')
                obj.SOS_b=obj.F.Numerator;
                obj.G_a=1;
            elseif isprop(obj.F,'sosMatrix')
                obj.SOS_b=obj.F.sosMatrix;
                obj.G_a=obj.F.scaleValues;
            else
                disp('Created filter should contain either a Numerator on SOS matrix, try again');
                [obj]=GUIdesign(obj);
            end
            obj.currentDesignedFilter=obj.dF.Response;
            obj.filteredSamplingFrequency=obj.samplingFrequency;
        end
    end
    methods (Hidden)
        function [obj]=filterData(samplingFrequency)
            %create filter object
            %Usage: [fil]=getFilteredData(samplingFrequency);
            %Input : samplingFrequency - the sampling frequency of the data
            %Output: fil - the filter obj
            %Assuming data format [ch x trials x samples]
            if nargin==1
                obj.samplingFrequency=samplingFrequency(1);
            elseif nargin==0
                obj.samplingFrequency=1;
                disp('Filter initialized with sampling freq. of 1Hz');
            end
        end
    end
end