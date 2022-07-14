function varargout=NlxCSC2mat(fname,FieldSelectionFlags,HeaderExtractionFlag,ExtractionMode,ExtractionModeVector)
% see help for Nlx2MatCSC. All functionalities are the same as in Nlx2MatCSC
% except for ExtractionMode = 3 or 5 which are not supported.

Fgetoutput=logical(FieldSelectionFlags);
outputnum=zeros(1,numel(Fgetoutput)+1);
for o=1:numel(Fgetoutput)
    outputnum(o)=sum(Fgetoutput(1:o));
end
outputnum(end)=sum(Fgetoutput)+1;

fid=fopen(fname,'r');
if fid==-1
    errordlg('file not found');
end



if sum(Fgetoutput)>0
    % go to end of file to get filesize.
    % note, that this info is also in the header... maybe in future do a double
    % check of the header info and the file size.
    fseek(fid, 0, 1);
    pos=ftell(fid);
    % the file is a 16K header with a bunch of CSC records afterwards.
    % each CSC record is 64+32+32+32+512*16 bits = 1044 bytes
    total_recs=(pos-16384)/1044;

    if mod(total_recs,1)>0
        warning('some bad records')
    end
    fseek(fid, 16384, 'bof'); %skip header

    switch ExtractionMode
        case 1
            start_rec=0;
            end_rec=total_recs;
        case 2
            start_rec=min(total_recs,ExtractionModeVector(1));
            end_rec=min(total_recs,ExtractionModeVector(2));
        case 3
            errordlg('ExtractionMode not supported : use Nlx2MatCSC instead');
        case 4
            fseek(fid, 16384, 'bof');
            ts0=fread(fid, 1, 'int64');
            channel0=fread(fid, 1, 'int32');
            sampFreq0=fread(fid, 1, 'int32');
            fseek(fid, 16384, 'bof');

            recstart=floor((ExtractionModeVector(1)-ts0)*sampFreq0*1e-6/512)+1;
            recend=floor((ExtractionModeVector(2)-ts0)*sampFreq0*1e-6/512)+1;
            start_rec=max(1,recstart);
            end_rec=max(1,recend);
            end_rec=min(end_rec,total_recs);
        case 5
            fseek(fid, 16384, 'bof');
            ts0=fread(fid, total_recs, 'int64',1044-8);
            ch0=fread(fid, 1, 'int32');
            Fs=fread(fid, 1000, 'int32',1044-4);
            
            errordlg('ExtractionMode not supported : use Nlx2MatCSC instead');        
    end

    nbrecs=end_rec-start_rec+1;

    ts=zeros(nbrecs, 1);
    channel=ts;
    sampFreq=ts;
    nValSamp=ts;
    data=zeros(nbrecs, 512);
    
    fseek(fid, 16384+1044*(start_rec-1), 'bof');
    for recX=1:nbrecs
        ts(recX)=fread(fid, 1, 'int64');
        channel(recX)=fread(fid, 1, 'int32');
        sampFreq(recX)=fread(fid, 1, 'int32');
        nValSamp(recX)=fread(fid, 1, 'int32');
        data(recX,:)=fread(fid,512,'int16');
    end
    fclose(fid);

    if Fgetoutput(1)
        varargout{1}=ts; 
    end
    if Fgetoutput(2)
        varargout{outputnum(2)}=channel; 
    end
    if Fgetoutput(3)
        varargout{outputnum(3)}=sampFreq; 
    end
    if Fgetoutput(4)
        varargout{outputnum(4)}=NumberOfValidSamples; 
    end
    if Fgetoutput(5)
        varargout{outputnum(5)}=data; 
    end

    if ~numel(unique(channel))==1
        warning('wierd channel information')        
    end

    if ~numel(unique(sampFreq))==1
        warning('wierd sampling frequency infomation')        
    end
end
end