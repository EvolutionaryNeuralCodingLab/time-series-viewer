load('layout_40_16x2_FlexLin.mat');
padSize=[12 35];

Enp(:,24)=[199;130];
Enp(:,25)=[205;170];
Enp(:,27)=[211;210];
Enp(:,26)=[217;250];
Enp(:,22)=[223;290];
Enp(:,23)=[229;330];
Enp(:,28)=[235;370];
Enp(:,29)=[241;410];
Enp(:,20)=[247;450];
Enp(:,21)=[253;490];
Enp(:,31)=[259;530];
Enp(:,30)=[265;570];
Enp(:,18)=[271;610];
Enp(:,19)=[277;650];
Enp(:,1)=[283;690];
Enp(:,32)=[289;730];

Enp(:,16)=[306;745];
Enp(:,17)=[311;705];
Enp(:,14)=[317;665];
Enp(:,15)=[323;625];
Enp(:,3)=[329;585];
Enp(:,2)=[335;545];
Enp(:,12)=[341;505];
Enp(:,13)=[347;465];
Enp(:,4)=[353;425];
Enp(:,5)=[359;385];
Enp(:,10)=[365;345];
Enp(:,11)=[371;305];
Enp(:,7)=[377;265];
Enp(:,6)=[383;225];
Enp(:,8)=[389;185];
Enp(:,9)=[395;145];

save('layout_40_16x2_FlexLin.mat','En','Ena','Enp');

%%


hight=600;
width=105;
dx=width/15;
dy=hight/15;

linearOrder=[24 25 27 26 22 23 28 29 20 21 31 30 18 19 1 32 16 17 14 15 3 2 12 13 4 5 10 11 7 6 8 9];

tmpPos=[0;0];
for i=1:16
    pElec=find(linearOrder==i);
    Enp(:,pElec)=tmpPos+[dx;-dy];
    tmpPos=Enp(:,pElec);
end

tmpPos(1)=tmpPos(1)+15-dx;
tmpPos(2)=tmpPos(2)-dy;
for i=17:32
    pElec=find(linearOrder==i);
    Enp(:,pElec)=tmpPos+[dx;dy];
    tmpPos=Enp(:,pElec);
end
Enp=Enp-min(Enp')'*ones(1,32);

save('layout_40_16x2_FlexLin.mat','En','Ena','Enp');


