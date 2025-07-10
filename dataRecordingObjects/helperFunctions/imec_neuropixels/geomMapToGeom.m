% =========================================================
% Parse snsGeomMap for XY coordinates
%
    function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta)

        C = textscan(meta.snsGeomMap, '(%d:%d:%d:%d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        shankInd = double(cell2mat(C(1)));
        xCoord = double(cell2mat(C(2)));
        yCoord = double(cell2mat(C(3)));
        connected = double(cell2mat(C(4)));

        % parse header for number of shanks
        geomStr = meta.snsGeomMap;
        headStr = extractBefore(geomStr,')(');
        headParts = split(headStr,',');
        nShank = str2double(headParts{2});
        shankWidth = str2double(headParts{4});
        shankPitch = str2double(headParts{3});
    end % geomMapToGeom