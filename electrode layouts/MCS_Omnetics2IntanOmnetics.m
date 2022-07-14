%% old
MCS=[1:2:31 2:2:32]; %IN
INTAN=[23:-1:8 24:31 0:7]+1;
IT2MCS(MCS)=INTAN;

%FlexLin32Ch
24 25 27 26 22 23 28 29 20 21 31 30 18 19 1 32 16 17 14 15 3 2 12 13 4 5 10 11 7 6 8 9

%FlexLin18Ch
24 25 23 26 21 29 19 30 1 17 3 14 4 12 7 10 9 8

%% Neu conversion
%electrode locations
En=[NaN   NaN   NaN    17    32   NaN   NaN   NaN;
   NaN   NaN   NaN    16     1   NaN   NaN   NaN;
   NaN   NaN   NaN    30     3   NaN   NaN   NaN;
   NaN   NaN    19    31     2    14   NaN   NaN;
   NaN   NaN    18   NaN   NaN    15   NaN   NaN;
   NaN   NaN    28   NaN   NaN     5   NaN   NaN;
   NaN    20    29   NaN   NaN     4    13   NaN;
   NaN    21   NaN   NaN   NaN   NaN    12   NaN;
   NaN    26   NaN   NaN   NaN   NaN     7   NaN;
    23    27   NaN   NaN   NaN   NaN     6    10;
    22   NaN   NaN   NaN   NaN   NaN   NaN    11;
    24   NaN   NaN   NaN   NaN   NaN   NaN     9;
    25   NaN   NaN   NaN   NaN   NaN   NaN     8];

Enp(:,25)=[199;130];
Enp(:,24)=[205;170];
Enp(:,22)=[211;210];
Enp(:,23)=[217;250];
Enp(:,27)=[223;290];
Enp(:,26)=[229;330];
Enp(:,21)=[235;370];
Enp(:,20)=[241;410];
Enp(:,29)=[247;450];
Enp(:,28)=[253;490];
Enp(:,18)=[259;530];
Enp(:,19)=[265;570];
Enp(:,31)=[271;610];
Enp(:,30)=[277;650];
Enp(:,16)=[283;690];
Enp(:,17)=[289;730];

Enp(:,32)=[306;745];
Enp(:,1)=[311;705];
Enp(:,3)=[317;665];
Enp(:,2)=[323;625];
Enp(:,14)=[329;585];
Enp(:,15)=[335;545];
Enp(:,5)=[341;505];
Enp(:,4)=[347;465];
Enp(:,13)=[353;425];
Enp(:,12)=[359;385];
Enp(:,7)=[365;345];
Enp(:,6)=[371;305];
Enp(:,10)=[377;265];
Enp(:,11)=[383;225];
Enp(:,9)=[389;185];
Enp(:,8)=[395;145]; %high end- hig

Enp(2,:)=800-Enp(2,:);
Enp(1,:)=Enp(1,:)-190;

Ena=cellfun(@(x) num2str(x),mat2cell(En,ones(1,size(En,1)),ones(1,size(En,2))),'UniformOutput',0);
save layout_40_16x2_FlexLin Enp En Ena;

