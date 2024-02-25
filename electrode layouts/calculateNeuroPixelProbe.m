En %chLayoutNumbers %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid
Ena %chLayoutNames %(Cell MxN)The layout of the channel names in physical space arranged in an M by N grid
electrodePitch %electrodePitch % distance between electrodes (not critical)
Enp %chLayoutPositions % (1xN or 2xN or 3xN) array of electrode position in [x or x,y or x,y,z]
layoutName %layoutName

load neuropixPhase3B1_kilosortChanMap.mat; %loaded from the configFiles folder in kilosort
layoutName='NP3B1';
Enp=[xcoords';ycoords'];
electrodePitch=[16 20];
x=(xcoords-min(xcoords))/electrodePitch(1)+1;
y=(ycoords-min(ycoords))/electrodePitch(2)+1;
En=nan(max(y),max(x));
En(sub2ind(size(En),y,x))=chanMap;
Ena=cellfun(@(x) num2str(x),mat2cell(En-1,ones(1,size(En,1)),ones(1,size(En,2))),'UniformOutput',0);
save layout_20_NP3B1 Ena En Enp electrodePitch layoutName;

load neuropixPhase3B2_kilosortChanMap.mat; %loaded from the configFiles folder in kilosort
layoutName='NP3B2';
Enp=[xcoords';ycoords'];
electrodePitch=[32 20];
x=(xcoords-min(xcoords))/electrodePitch(1)+1;
y=(ycoords-min(ycoords))/electrodePitch(2)+1;
En=nan(max(y),max(x));
En(sub2ind(size(En),y,x))=chanMap;
Ena=cellfun(@(x) num2str(x),mat2cell(En-1,ones(1,size(En,1)),ones(1,size(En,2))),'UniformOutput',0);
save layout_20_NP3B2 Ena En Enp electrodePitch layoutName;
