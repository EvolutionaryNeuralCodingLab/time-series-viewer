%%% CatGT command line

%Important To take into account before running:
%1. Remember to separate recordings into a new 'Insertionx' folder where x is
%   the number of the insertion.
%2. The name of the folder holding all insertions must be the same name
%   given to all acquisition events in spike GLX (each insertion is one
%   acquisition event, do not change names between acquisitions in an experimental day). 

%example variables:

% base_dir = "\\132.66.45.127\data\Simon\Anesthesia_exp\Experiment_29_09_22";
% run_dir = "\\132.66.45.127\data\Simon\SpikeGLX\"; %Folder that has catgt and tprime subfolders
% insertion = "4";


%% CatGT function
function CatGT_Tprime(base_dir, run_dir, insertion)
    
    %0.0 Set directory to catGT directory
    cd(run_dir+"CatGT-win")

    %0.1 Find experiment file name from base_dir
    out=regexp(base_dir,'\','split');
    exp = string(out(end));
    
    %0.2 Find number of stimulus in insertion file
    file = dir (base_dir + "\Insertion" + insertion);
    filenames = {file.name};
    num = string(sum( ~cellfun(@isempty, strfind(filenames, 'Experiment')))-1);

    %1. Create command line filter + concatenation ap
    cmndAP = "CatGT -dir="+ base_dir + "\Insertion" + insertion + " -run=" + exp +"_" + insertion ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -ap" + " -prb=0 -zerofillmax=0"...
        + " -apfilter=butter,12,300,9000 -gfix=0.40,0.10,0.02 -loccar=4,32 -dest="...
        + base_dir + "\Insertion" + insertion;
    
    status = 1;
    
    %1.1. execute command
    status = system(cmndAP);
   
    while status == 1
    %waiting for command to complete
    end

    disp("AP extraction and concatenation completed")
    
    %2.0. Create lf folder
    path = base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0";
    
    cd(path)
    
    if ~exist('sync_events','dir')
        mkdir sync_events
    end
    
    %Go to TPrime dir
    cd(run_dir + "TPrime-win")

    %2. Create command line filter + concatenation lf
    cmndLF = "CatGT -dir="+ base_dir + "\Insertion" + insertion + " -run=" + exp +"_" + insertion ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -lf" + " -prb=0 -zerofillmax=0"...
        + " -dest="...
        + base_dir + "\Insertion" + insertion;

    status = 1;
    
    %2.1. execute command
    status = system(cmndLF);
   
    while status == 1
    %waiting for command to complete
    end
    
    disp("LF extraction and concatenation completed")

    %3. Create command line for concatenation nidaq
    cmndNIc = "CatGT -dir="+ base_dir + "\Insertion" + insertion + " -run=" + exp +"_" + insertion ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -ni" + " -prb=0 -zerofillmax=0"...
        + " -dest="...
        + base_dir + "\Insertion" + insertion;

    status = 1;
    
    %3.1. execute command
    status = system(cmndNIc);
   
    while status == 1
    %waiting for command to complete
    end
    
    disp("NI concatenation completed")

    %4. Create command line for extracting nidaq
    cmndNIe = "CatGT -dir="+ base_dir + "\Insertion" + insertion + " -run=catgt_" + exp +"_" + insertion ...
        + " -g=0 -ni -prb=0 -t=cat -no_tshift "...
        + "-xd=0,0,-1,0,0 "...
        + "-xd=0,0,-1,1,0 "... 
        + "-xd=0,0,-1,2,0 "...
        + "-xd=0,0,-1,3,0 "...
        + "-xid=0,0,-1,1,0 "...
        + "-xid=0,0,-1,2,0 "... 
        + "-dest="...
        + base_dir + "\Insertion" + insertion;

    status = 1;
    
    %4.1. execute command
    status = system(cmndNIe);
   
    while status == 1
    %waiting for command to complete
    end
    
    disp("NI extraction completed")
    
    %5.0. Create sync folder
    path = base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0";
    
    cd(path)
    
    if ~exist('sync_events','dir')
        mkdir sync_events
    end
    
    %Go to TPrime dir
    cd(run_dir + "TPrime-win")

    %5.1. Create TPrime command
    cmndTPrime = "TPrime -syncperiod=1.0 -tostream=" + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.imec0.ap.xd_384_6_500.txt "...
    + "-fromstream=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_8_0_500.txt "...
    + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_8_1_0.txt,"...
    + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_1.txt "...
    + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_8_2_0.txt,"...
    + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_2.txt "...
    + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_8_3_0.txt,"...
    + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_3.txt "...
    + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xid_8_1_0.txt,"...
    + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_1inv.txt "...
    + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xid_8_2_0.txt,"...
    + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_2inv.txt";
    
    status = 1;
    
    %5.2. execute command
    status = system(cmndTPrime);
   
    while status == 1
    %waiting for command to complete
    end
    
    disp("Event syncing completed, look for files in sync_events folder")
  
end




