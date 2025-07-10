% =========================================================
% Return geometry paramters for supported probe types
% These are used to calculate positions from metadata
% that includes only ~snsShankMap
%
    function geom = getGeomParams(meta)
        % create map
        geomTypeMap = makeTypeMap();

        % get probe part number; if absent, this is a 3A
        if isfield(meta,'imDatPrb_pn')
            pn = meta.imDatPrb_pn;
        else
            pn = '3A';
        end

        if geomTypeMap.isKey(pn)
            geomType = geomTypeMap(pn);
        else
            fprintf('unsupported probe part number\n');
            return;
        end

        switch geomType
            case 'np1_stag_70um'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 11;
                geom.horzPitch = 32;
                geom.vertPitch = 20;
                geom.rowsPerShank = 480;
                geom.elecPerShank = 960;
            case 'nhp_lin_70um'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 27;
                geom.horzPitch = 32;
                geom.vertPitch = 20;
                geom.rowsPerShank = 480;
                geom.elecPerShank = 960;
            case 'nhp_stag_125um_med'
                geom.nShank = 1;
                geom.shankWidth = 125;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 11;
                geom.horzPitch = 87;
                geom.vertPitch = 20;
                geom.rowsPerShank = 1368;
                geom.elecPerShank = 2496;
            case 'nhp_stag_125um_long'
                geom.nShank = 1;
                geom.shankWidth = 125;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 11;
                geom.horzPitch = 87;
                geom.vertPitch = 20;
                geom.rowsPerShank = 2208;
                geom.elecPerShank = 4416;
            case 'nhp_lin_125um_med'
                geom.nShank = 1;
                geom.shankWidth = 125;
                geom.shankPitch = 0;
                geom.even_xOff = 11;
                geom.odd_xOff = 11;
                geom.horzPitch = 103;
                geom.vertPitch = 20;
                geom.rowsPerShank = 1368;
                geom.elecPerShank = 2496;
            case 'nhp_lin_125um_long'
                geom.nShank = 1;
                geom.shankWidth = 125;
                geom.shankPitch = 0;
                geom.even_xOff = 11;
                geom.odd_xOff = 11;
                geom.horzPitch = 103;
                geom.vertPitch = 20;
                geom.rowsPerShank = 2208;
                geom.elecPerShank = 4416;
            case 'uhd_8col_1bank'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 14;
                geom.odd_xOff = 14;
                geom.horzPitch = 6;
                geom.vertPitch = 6;
                geom.rowsPerShank = 48;
                geom.elecPerShank = 384;
            case 'uhd_8col_16bank'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 14;
                geom.odd_xOff = 14;
                geom.horzPitch = 6;
                geom.vertPitch = 6;
                geom.rowsPerShank = 768;
                geom.elecPerShank = 6144;
            case 'np2_ss'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 27;
                geom.horzPitch = 32;
                geom.vertPitch = 15;
                geom.rowsPerShank = 640;
                geom.elecPerShank = 1280;
            case 'np2_4s'
                geom.nShank = 4;
                geom.shankWidth = 70;
                geom.shankPitch = 250;
                geom.even_xOff = 27;
                geom.odd_xOff = 27;
                geom.horzPitch = 32;
                geom.vertPitch = 15;
                geom.rowsPerShank = 640;
                geom.elecPerShank = 1280;
            case 'NP1120'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 6.75;
                geom.odd_xOff = 6.75;
                geom.horzPitch = 4.5;
                geom.vertPitch = 4.5;
                geom.rowsPerShank = 192;
                geom.elecPerShank = 384;
            case 'NP1121'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 6.25;
                geom.odd_xOff = 6.25;
                geom.horzPitch = 3;
                geom.vertPitch = 3;
                geom.rowsPerShank = 384;
                geom.elecPerShank = 384;
            case 'NP1122'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 12.5;
                geom.odd_xOff = 12.5;
                geom.horzPitch = 3;
                geom.vertPitch = 3;
                geom.rowsPerShank = 24;
                geom.elecPerShank = 384;
            case 'NP1123'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 10.25;
                geom.odd_xOff = 10.25;
                geom.horzPitch = 4.5;
                geom.vertPitch = 4.5;
                geom.rowsPerShank = 32;
                geom.elecPerShank = 384;
            case 'NP1300'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 11;
                geom.odd_xOff = 11;
                geom.horzPitch = 48;
                geom.vertPitch = 20;
                geom.rowsPerShank = 480;
                geom.elecPerShank = 960;
            case 'NP1200'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 11;
                geom.horzPitch = 32;
                geom.vertPitch = 20;
                geom.rowsPerShank = 64;
                geom.elecPerShank = 128;
            case 'NXT3000'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 53;
                geom.odd_xOff = 53;
                geom.horzPitch = 0;
                geom.vertPitch = 15;
                geom.rowsPerShank = 128;
                geom.elecPerShank = 128;
            otherwise
                % shouldn't see this case
                fprintf('unsupported probe part number\n');
                return;
        end
    end %end  getGeomParam