load neuronexuxIntanConversion32;

dX=200;
dY=50;

probe=[(1:8)' (9:16)' (17:24)' (25:32)'];
probeX=ones(8,1)*(0:dX:600);
probeY=flipud((0:dY:350)'*ones(1,4));

probe=flipud(probe);
probeX=flipud(probeX);
probeY=flipud(probeY);

probeIntan=NeuroNexus2Intan(probe);
probeIntanFlipped=NeuroNexus2IntanFlipped(probe);

probeName='NeuroNexus: A8x4-2mm-50-200-177-A32';
En=probeIntan;
Ena=cellfun(@(x) num2str(x),mat2cell(En,ones(1,size(En,1)),ones(1,size(En,2))),'UniformOutput',0);
for i=unique(probeIntan)'
    p=find(probe==i);
    Enp(:,i)=[probeX(p);probeY(p)];
end
save layout_50_8x4_EdgeProbe probeName En Ena Enp;

probeName='NeuroNexus: A8x4-2mm-50-200-177-A32 - Flipped';
En=probeIntanFlipped;
Ena=cellfun(@(x) num2str(x),mat2cell(En,ones(1,size(En,1)),ones(1,size(En,2))),'UniformOutput',0);
for i=unique(probeIntanFlipped)'
    p=find(probe==i);
    Enp(:,i)=[probeX(p);probeY(p)];
end
save layout_50_8x4_EdgeProbeFlipped probeName En Ena Enp;
