 function M = makeTypeMap()
        % many part numbers have the same geometry parameters ;
        % make a map that pairs geometry type (value) with probe part number (key)
        M = containers.Map('KeyType','char','ValueType','char');

        M('3A') = 'np1_stag_70um';
        M('PRB_1_4_0480_1') = 'np1_stag_70um';
        M('PRB_1_4_0480_1_C') = 'np1_stag_70um';
        M('NP1010') = 'np1_stag_70um';
        M('NP1011') = 'np1_stag_70um';
        M('NP1012') = 'np1_stag_70um';
        M('NP1013') = 'np1_stag_70um';

        M('NP1015') = 'nhp_lin_70um';
        M('NP1015') = 'nhp_lin_70um';
        M('NP1016') = 'nhp_lin_70um';
        M('NP1017') = 'nhp_lin_70um';

        M('NP1020') = 'nhp_stag_125um_med';
        M('NP1021') = 'nhp_stag_125um_med';
        M('NP1030') = 'nhp_stag_125um_long';
        M('NP1031') = 'nhp_stag_125um_long';

        M('NP1022') = 'nhp_lin_125um_med';
        M('NP1032') = 'nhp_lin_125um_long';

        M('NP1100') = 'uhd_8col_1bank';
        M('NP1110') = 'uhd_8col_16bank';

        M('PRB2_1_2_0640_0') = 'np2_ss';
        M('PRB2_1_4_0480_1') = 'np2_ss';
        M('NP2000') = 'np2_ss';
        M('NP2003') = 'np2_ss';
        M('NP2004') = 'np2_ss';

        M('PRB2_4_2_0640_0') = 'np2_4s';
        M('PRB2_4_4_0480_1') = 'np2_4s';
        M('NP2010') = 'np2_4s';
        M('NP2013') = 'np2_4s';
        M('NP2014') = 'np2_4s';

        M('NP1120') = 'NP1120';
        M('NP1121') = 'NP1121';
        M('NP1122') = 'NP1122';
        M('NP1123') = 'NP1123';
        M('NP1300') = 'NP1300';

        M('NP1200') = 'NP1200';
        M('NXT3000') = 'NXT3000';
    end % makeTypeMap