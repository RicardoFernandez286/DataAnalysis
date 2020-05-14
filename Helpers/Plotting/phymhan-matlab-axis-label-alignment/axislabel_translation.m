function axislabel_translation(h_a,trans_pnt,trans_vec,trans_units)
% This function moves axis labels to a proper distance from the axes. It is
% called by function 'align_axislabels', which first, rotates the x, y and 
% z labels.
% If you just want to move the labels but not rotate them, you have to 
% specify the position where to put the labels and the direction along 
% which you move labels away from the axes.
%
% Input arguments: 
%  trans_units: 'pixels', move labels away from axes in pixels;
%  'characters', move labels away from axes in characters;
%  trans_vec: translation vector, a 2-by-3 matrix, each column defines the 
%  direction along which you move labels away from the axes, the vector 
%  should be normalized;
%  trans_pnt: translation point, a 3-by-3 matrix, each column defines the 
%  position (in the original 3D space rather than in the 2D canvas space) 
%  where you put the axis labels;
%  h_a: handle of the axes.
%
% Apr/02/2015

global AXISALIGN_TRANS_A AXISALIGN_TRANS_B
if isempty(AXISALIGN_TRANS_A), AXISALIGN_TRANS_A = 1; end
if isempty(AXISALIGN_TRANS_B), AXISALIGN_TRANS_B = 1; end

if nargin < 4
    trans_units = 'pixels';
end
vers = version();
vers = str2double(vers(1:3));
h_xlb = get(h_a,'xlabel');
h_ylb = get(h_a,'ylabel');
h_zlb = get(h_a,'zlabel');
if strcmpi(trans_units,'pixels')
    trans_units = 1;
    % Modify the method by which you move labels here
    %======================================%
    char_pixel_len = get(h_xlb,'FontSize')*0.55; % 0.67
    base_pixel_len = get(h_xlb,'FontSize')*2.0; % 2.0
    %======================================%
elseif strcmpi(trans_units,'character')
    trans_units = 2;
    % Modify the method by which you move labels here
    %======================================%
    % pos = pos+trans_vec'.*aniso_ratio*(A*char_len+B)
    A = AXISALIGN_TRANS_A;
    B = AXISALIGN_TRANS_B;
    aniso_ratio = [2, 1];
    %======================================%
end
%% Move XLabel
set(h_xlb,'HorizontalAlignment','center','VerticalAlignment','Middle','Units','data');
set(h_xlb,'Position',trans_pnt(:,1))
char_len_x = label_max_len('x');
if trans_units == 2
    set(h_xlb,'units','characters');
    pos_xlb = get(h_xlb,'Position');
    pos_xlb(1:2) = pos_xlb(1:2)+trans_vec(:,1)'.*aniso_ratio*(A*char_len_x+B);
else
    set(h_xlb,'Units','Pixels');
    pos_xlb = get(h_xlb,'Position');
    pos_xlb(1:2) = pos_xlb(1:2)+trans_vec(:,1)'*(base_pixel_len+char_len_x*char_pixel_len);
end
set(h_xlb,'Position',pos_xlb)
%% Move YLabel
set(h_ylb,'HorizontalAlignment','center','VerticalAlignment','Middle','Units','data');
set(h_ylb,'Position',trans_pnt(:,2))
char_len_y = label_max_len('y');
if trans_units == 2
    set(h_ylb,'units','characters');
    pos_ylb = get(h_ylb,'Position');
    pos_ylb(1:2) = pos_ylb(1:2)+trans_vec(:,2)'.*aniso_ratio*(A*char_len_y+B);
else
    set(h_ylb,'Units','Pixels');
    pos_ylb = get(h_ylb,'Position');
    pos_ylb(1:2) = pos_ylb(1:2)+trans_vec(:,2)'*(base_pixel_len+char_len_y*char_pixel_len);
end
set(h_ylb,'Position',pos_ylb)
%% Move ZLabel
set(h_zlb,'HorizontalAlignment','center','VerticalAlignment','Middle','Units','data');
set(h_zlb,'Position',trans_pnt(:,3))
char_len_z = label_max_len('z');
if trans_units == 2
    set(h_zlb,'units','characters');
    pos_zlb = get(h_zlb,'Position');
    pos_zlb(1:2) = pos_zlb(1:2)+trans_vec(:,3)'.*aniso_ratio*(A*char_len_z+B);
else
    set(h_zlb,'Units','Pixels');
    pos_zlb = get(h_zlb,'Position');
    pos_zlb(1:2) = pos_zlb(1:2)+trans_vec(:,3)'*(base_pixel_len+char_len_z*char_pixel_len);
end
set(h_zlb,'Position',pos_zlb)

    function maxlen = label_max_len(ax)
        ticklabel = get(h_a, [ax, 'TickLabel']);
        if vers < 8.4
            maxlen = size(ticklabel,2);
        else
            maxlen = 0;
            for k_f = 1:length(ticklabel)
                if length(ticklabel{k_f}) > maxlen
                    maxlen = length(ticklabel{k_f});
                end
            end
        end
    end
end
