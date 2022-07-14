function [ch , t]= parseHeka(data,ser,chNumbers)
% [ch , t]= parseHeka(data,ser,chNumbers)
% Function purpose : Parse data form patch clamp heka files in .mat format
%
% Function recives :    data - a structure in matlab format
%                       ser - the series to extract
%                       chNumbers - the channels to extract (in case more than one channel was recorded)
%
% Function give back :  ch - activity traces
%                       t - time vector

%
% Last updated : 24/10/14
nCh=numel(chNumbers);

if numel(ser)==1
    if isfield(data.series,'trace')
        for n=1:nCh;
            ch(n).vm = data.series(ser).trace(1,chNumbers(n)).data;
            ch(n).mvm = mean(ch(n).vm,2); % mean Vm trace
            ch(n).units = data.series(ser).trace(1,chNumbers(n)).header.adc.Units;
        end
        t.dt =data.series(ser).trace(chNumbers(1)).timeInfo.dt;
        t.s = data.series(ser).trace(chNumbers(1)).timeInfo.dt:data.series(ser).trace(chNumbers(1)).timeInfo.dt:length(ch(n).vm(:,1))*data.series(ser).trace(chNumbers(1)).timeInfo.dt;
        t.ms = t.s.*1000;
        
    elseif isfield(data.series,'traces')
        for n=1:nCh;
            ch(n).vm = data.series(ser).traces(1,chNumbers(n)).values';
            ch(n).mvm = mean(ch(n).vm,2); % mean Vm trace
            ch(n).units = data.series(ser).traces(1,chNumbers(n)).unit;
        end
        
        t.dt =data.series(ser).traces(chNumbers(1)).timeVector; t.dt = t.dt(2);
        t.s = t.dt:t.dt:length(ch(n).vm(:,1))*t.dt;
        t.ms = t.s.*1000;
    end
end


if numel(ser)>1
    if isfield(data.series,'trace')
        for z=1:numel(ser)
            tmp_time(z) = length(data.series(ser(z)).trace(1,1).data(:,1)); %change to the first channel in the series as before chNumbers(1)
            tmp_trials(z) = length(data.series(ser(z)).trace(1,1).data(1,:));
        end
        ser_maxdatapoints= max(tmp_time); % index of series with the maximal time datapoints
        maxtrials= sum(tmp_trials); % sum of all the trials
        
        means = [];
        for z=1:numel(ser)
            for n=1:nCh;
                parse(z).ch(n).vm = data.series(ser(z)).trace(1,chNumbers(n)).data;
                parse(z).ch(n).mvm = mean(parse(z).ch(n).vm,2); % mean Vm trace
                parse(z).ch(n).trialmeans = [means,mean(parse(z).ch(n).vm,1)];
                parse(z).ch(n).units = data.series(ser(z)).trace(1,chNumbers(n)).header.adc.Units;
                %size(data.series(ser(ser_maxdatapoints)).trace(1,1).data) % size of the series with the maximal data point
            end
        end
        
        % create final big channel matrix
        tmp_mat = zeros(ser_maxdatapoints,maxtrials,nCh);
        tmp_meanmat = zeros(ser_maxdatapoints,numel(ser),nCh);
        cumsum_trials = [0 ,cumsum(tmp_trials)];
        tmp_units = cell(numel(ser),nCh);
        for i = 1:numel(ser)
            for n=1:nCh;
                tmp_mat(:,cumsum_trials(i)+1:cumsum_trials(i+1),n) = repmat(parse(i).ch(n).trialmeans,ser_maxdatapoints,1);
                tmp_meanmat(:,i,n) = repmat(mean(parse(i).ch(n).mvm),ser_maxdatapoints,1);
                tmp_mat(1:length(parse(i).ch(n).vm),cumsum_trials(i)+1:cumsum_trials(i+1),n) = parse(i).ch(n).vm;
                tmp_meanmat(1:length(parse(i).ch(n).mvm),i,n) = parse(i).ch(n).mvm;
                tmp_units(i,n)=parse(i).ch(n).units;
            end
        end
        
        for z=1:numel(ser)
            for n=1:nCh
                ch(n).vm = tmp_mat(:,:,n);
                ch(n).mvm = mean(ch(n).vm,2);
                ch(n).ind_mvm = tmp_meanmat(:,:,n);
                ch(n).units = tmp_units(1,n);
            end
        end
        
        t.dt =data.series(ser(1)).trace(chNumbers(1)).timeInfo.dt;
        t.s = t.dt:t.dt:length(tmp_mat(:,1,1))*t.dt;
        t.ms = t.s.*1000;
        
    elseif isfield(data.series,'traces')
        for z=1:numel(ser)
            tmp_time(z) = length(data.series(ser(z)).traces(1,1).values(1,:)'); %change to the first channel in the series as before chNumbers(1)
            tmp_trials(z) = length(data.series(ser(z)).traces(1,1).values(1,:)');
        end
        ser_maxdatapoints= max(tmp_time); % index of series with the maximal time datapoints
        maxtrials= sum(tmp_trials); % sum of all the trials
        
        means = [];
        for z=1:numel(ser)
            for n=1:nCh;
                parse(z).ch(n).vm = data.series(ser(z)).traces(1,chNumbers(n)).data;
                parse(z).ch(n).mvm = mean(parse(z).ch(n).vm,2); % mean Vm trace
                parse(z).ch(n).trialmeans = [means,mean(parse(z).ch(n).vm,1)];
                parse(z).ch(n).units = data.series(ser(z)).traces(1,chNumbers(n)).header.adc.Units;
                %size(data.series(ser(ser_maxdatapoints)).traces(1,1).data) % size of the series with the maximal data point
            end
        end
        
        % create final big channel matrix
        tmp_mat = zeros(ser_maxdatapoints,maxtrials,nCh);
        tmp_meanmat = zeros(ser_maxdatapoints,numel(ser),nCh);
        cumsum_trials = [0 ,cumsum(tmp_trials)];
        tmp_units = cell(numel(ser),nCh);
        for i = 1:numel(ser)
            for n=1:nCh;
                tmp_mat(:,cumsum_trials(i)+1:cumsum_trials(i+1),n) = repmat(parse(i).ch(n).trialmeans,ser_maxdatapoints,1);
                tmp_meanmat(:,i,n) = repmat(mean(parse(i).ch(n).mvm),ser_maxdatapoints,1);
                tmp_mat(1:length(parse(i).ch(n).vm),cumsum_trials(i)+1:cumsum_trials(i+1),n) = parse(i).ch(n).vm;
                tmp_meanmat(1:length(parse(i).ch(n).mvm),i,n) = parse(i).ch(n).mvm;
                tmp_units(i,n)=parse(i).ch(n).units;
            end
        end
        
        for z=1:numel(ser)
            for n=1:nCh
                ch(n).vm = tmp_mat(:,:,n);
                ch(n).mvm = mean(ch(n).vm,2);
                ch(n).ind_mvm = tmp_meanmat(:,:,n);
                ch(n).units = tmp_units(1,n);
            end
        end
        
        t.dt =data.series(ser(1)).traces(chNumbers(1)).timeInfo.dt;
        t.s = t.dt:t.dt:length(tmp_mat(:,1,1))*t.dt;
        t.ms = t.s.*1000;
    end
end