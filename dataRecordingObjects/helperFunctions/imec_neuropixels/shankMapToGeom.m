% =========================================================
% Get XY coordinates from snsShankMap plus hard coded geom values
%
    function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta)
        % get number of saved AP channels (some early metadata files have a
        % SYNC entry in the snsChanMap
        [nchan,~,~] = ChannelCountsIM(meta);

        C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        shankInd = double(cell2mat(C(1)));
        colInd = double(cell2mat(C(2)));
        rowInd = double(cell2mat(C(3)));
        connected = double(cell2mat(C(4)));

        % trim these to the number of saved channels
        shankInd = shankInd(1:nchan);
        colInd = colInd(1:nchan);
        rowInd = rowInd(1:nchan);
        connected = connected(1:nchan);

        geom = getGeomParams(meta);

        oddRows = logical(mod(rowInd,2));
        evenRows = ~oddRows;
        xCoord = colInd*geom.horzPitch;
        xCoord(evenRows) = xCoord(evenRows) + geom.even_xOff ;
        xCoord(oddRows) = xCoord(oddRows) + geom.odd_xOff;
        yCoord = rowInd*geom.vertPitch;

        nShank = geom.nShank;
        shankWidth = geom.shankWidth;
        shankPitch = geom.shankPitch;
    end % shankMapToGeom