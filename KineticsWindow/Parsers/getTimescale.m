function powS = getTimescale(timescalestring)

switch timescalestring
    case 's'
        powS = 1;
    case 'ms'
        powS = -3;
    case {'us','\mus'}
        powS = -6;
    case 'ns'
        powS = -9;
    case 'ps'
        powS = -12;
    case 'fs'
        powS = -15;
    otherwise
        powS = 1;
end