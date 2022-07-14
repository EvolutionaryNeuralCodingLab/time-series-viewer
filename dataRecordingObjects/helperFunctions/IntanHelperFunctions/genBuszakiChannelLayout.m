
probe=[ 5    13    21    29
     4    12    20    28
     6    14    22    30
     3    11    19    27
     7    15    23    31
     2    10    18    26
     8    16    24    32
     1     9    17    25];
 adaptorToIntan=[ 31
    27
    22
    18
    28
    23
    21
    26
    29
    24
    20
    25
    30
    19
    32
    17
     1
    16
     3
    14
     9
    10
     8
     2
     7
    15
    11
    12
     6
    13
     5
     4
];


% s=1:2:7;
% En=nan(8,7);
% for shank=1:4
%     for p=1:8
%      En(p,s(shank))=adaptorToIntan(probe(p,shank));
%     end
% end
% 
% 
% letters=[{'A'},{[]},{ 'B'},{[]}, {'C'}, {[]},{'D'}];
% 
% for shank=[1:2:7]
%     for p=1:8
%         Ena{p,shank}=[letters{shank} int2str(p)];
%     end
% end

s=1:4;
En=nan(8,4);
for shank=1:4
    for p=1:8
     En(p,s(shank))=adaptorToIntan(probe(p,shank));
    end
end


letters=[{'A'},{ 'B'}, {'C'},{'D'}];

for shank=[1:4]
    for p=1:8
        Ena{p,shank}=[letters{shank} int2str(p)];
    end
end

save('layout_32_BuszakiProbe.mat','En','Ena')
