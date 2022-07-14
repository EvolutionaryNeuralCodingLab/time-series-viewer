Intan=[8:31 0:7]+1;
IntanFlipped=[24:31 0:7 8:23]+1;

Adapter=[23:2:31 19 17 21 11 15 13 1:2:9 10:-2:2 14 16 12 22 18 20 32:-2:24];
%Adapter=[9:-2:1 13 15 11 21 17 19 31:-2:23 24:2:32 20 18 22 12 16 14 2:2:10];

Adapter2Intan(Adapter)=Intan;
Adapter2IntanFlipped(IntanFlipped)=Adapter;

neuroNexus=[12 10 8 6 2 3 5 7 9 11 16 15 14 13 4 1 17 18 19 20 24 26 21 23 22 25 27 29 28 31 30 32];
adapter=[10:-1:1 11:1:16 22:-1:17 23:1:32];

NeuroNexus2Adapter(neuroNexux)=adapter;

NeuroNexus2Intan=Adapter2Intan(NeuroNexus2Adapter);
NeuroNexus2IntanFlipped=Adapter2IntanFlipped(NeuroNexus2Adapter);

save neuronexuxIntanConversion32 Adapter2Intan Adapter2IntanFlipped neuroNexus adapter NeuroNexus2Intan NeuroNexus2IntanFlipped;