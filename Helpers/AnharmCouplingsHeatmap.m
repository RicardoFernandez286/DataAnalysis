Nmodes = 4;
ModeLabels = ["A'(2)" "A''" "A'(1)" "C{\equiv}C"];

molec = 'NMe2';

switch molec
    case 'NMe2'
freqs = [...
1919.52
1942.86
2047.64
2214.33
3825.60
3874.80
4087.11
4410.94
3855.66
3953.72
4133.74
3971.61
4157.09
4261.84
];

    case 'Ph'
freqs = [...
1922.71
1945.21
2049.34
2254.69
3831.83
3879.38
4090.51
4491.78
3861.38
3958.33
4177.40
3975.65
4199.89
4303.99
];
end

%% Sort the data into the one and two exciton blocks
NH2 = Nmodes + nchoosek(Nmodes,2);
ExpectedN = Nmodes + NH2;

if length(freqs) ~= ExpectedN
    disp('Wrong number of modes')
    return
end

fund    = freqs(1:Nmodes);
ovcom   = freqs(Nmodes+1:end);

[ii,jj] = ndgrid(1:Nmodes);

H1 = fund*ones(Nmodes,1)'+(fund*ones(Nmodes,1)')';
H2 = zeros(Nmodes);

H2 = diag(ovcom(1:Nmodes),0);
H2(ii>jj) = ovcom(Nmodes+1:end);
H2 = triu(H2.',1) + tril(H2);

anh = H1-H2;



%% Plot and beautify
fh1 = figure(1);
clf(fh1);

anh(ii<jj) = nan; % To remove the upper block
% anh(ii==jj) = nan; % To remove diagonal

NanCol  = 0.9*[1 1 1];
LineCol = 0.9*[1 1 1];
h_map = heatmap(fh1,anh,'MissingDataLabel','','CellLabelFormat','%0.1f','MissingDataColor',NanCol,'FontSize',16);
colormap(othercolor('Blues9',Nmodes));
% colormap(othercolor('Reds9'));
% colormap(othercolor('Oranges9'));


warning('off','MATLAB:structOnObject');

hHeatmap = struct(h_map).Heatmap;
hGrid = struct(hHeatmap).Grid;
hGrid.ColorData = uint8(256*[LineCol 0.5]');

h_map.XData = ModeLabels;
h_map.YData = ModeLabels;