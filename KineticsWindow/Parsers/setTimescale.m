function TS = setTimescale(powS)
% This function outputs the correct timescale string for the appropriate power of ten from the seconds (powS)

switch powS
    case 0
        TS = 's';
    case -3
        TS = 'ms';
    case -6
        TS = '\mus';
    case -9
        TS = 'ns';
    case -12
        TS = 'ps';
    case -15
        TS = 'fs';
    otherwise
        if powS >= 0
            TS = 's';
        end
end

