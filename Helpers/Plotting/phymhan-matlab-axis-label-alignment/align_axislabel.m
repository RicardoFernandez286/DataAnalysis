function align_axislabel(~,ax)
% This function first rotates x, y and z labels to the direction of their 
% corresponding axes when a rotate operation has finished. So that labels 
% will be parallel to axes. Then, moves the axis labels to a proper 
% distance from the axes by calling 'axislabel_translation'.
% The function is used as an 'ActionPostCallback', which is a callback 
% function of a rotate3d object.
% It still works when the projection mode is perspective or when the data 
% aspect ratio is not [1 1 1].
%
% Example 1: % set ActionPostCallback to align_axislabel
%  h = rotate3d;
%  set(h,'ActionPostCallback',@align_axislabel);
%
% Example 2: % set WindowButtonMotionFcn to align_axislabel
%  h = rotate3d;
%  set(h,'ActionPreCallback',...
%      'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%  set(h,'ActionPostCallback',...
%      'set(gcf,''windowbuttonmotionfcn'','''')')
%
% Jul/13/2015, Adds a trans_mode option
% Feb/06/2015
% Apr/10/2016

% trans_mode controls the way we define trans_vec_x and trans_vec_y, which
% is an input argument passed to function axislabel_translation.
% Take x-axis for example. If trans_mode == 2, the translation direction is
% parallel to y-axis, and if trans_moede == 1, the translation direction is
% perpendicular to x-axis.
trans_mode = 1;

if isa(ax,'struct')
    h_a = ax.Axes;
elseif isa(ax,'matlab.graphics.axis.Axes')
    h_a = ax;
else
    h_a = gca;
end

proj_mode = get(h_a,'projection');
if strcmpi(proj_mode,'orthographic')
    proj_mode = 1;
else
    proj_mode = 2;
end
data_aspect_ratio = get(h_a,'dataaspectratio');
data_aspect_ratio = data_aspect_ratio(:);

set(h_a, 'DataAspectRatio', data_aspect_ratio); % optional, you can comment this line out

x_lim = get(h_a,'xlim');
y_lim = get(h_a,'ylim');
z_lim = get(h_a,'zlim');
cam_pos = get(h_a,'cameraposition');
cam_tar = get(h_a,'cameratarget');
cam_vec = cam_pos-cam_tar;
N = norm(cam_vec,2); % Near field
[az,el] = view; % get current azimuth and elevation angles
if abs(el) == 90
    az = atan2d(sind(az)*data_aspect_ratio(2),cosd(az)*data_aspect_ratio(1));
end % Adjust az when el == 90
az = az/180*pi;
el = pi/2-el/180*pi;
% Rotation Matrices
R_az = [cos(-az),-sin(-az),0,0;sin(-az),cos(-az),0,0;0,0,1,0;0,0,0,1];
R_el = [1,0,0,0;0,cos(-el),-sin(-el),0;0,sin(-el),cos(-el),0;0,0,0,1];
% Translation Matrices
T_tar = [1,0,0,-cam_tar(1);0,1,0,-cam_tar(2);0,0,1,-cam_tar(3);0,0,0,1];
T_cam = [1,0,0,0;0,1,0,0;0,0,1,-N;0,0,0,1];
pnt = [x_lim(1),y_lim(1),z_lim(1),1;...
    x_lim(2),y_lim(1),z_lim(1),1;...
    x_lim(2),y_lim(2),z_lim(1),1;...
    x_lim(1),y_lim(2),z_lim(1),1;...
    x_lim(1),y_lim(1),z_lim(2),1;...
    x_lim(2),y_lim(1),z_lim(2),1;...
    x_lim(2),y_lim(2),z_lim(2),1;...
    x_lim(1),y_lim(2),z_lim(2),1]';
pnt_cam = zeros(4,8);
if proj_mode == 1
    F = N+100; % Far field, the value would not influence the projection
    a = -(F+N)/(F-N);
    b = -2*F*N/(F-N);
    M_orth = [N,0,0,0;0,N,0,0;0,0,a,b;0,0,-1,0];
    pnt_cam = zeros(4,4);
    for kk = 1:8
        p = pnt(:,kk);
        p = T_tar*p;
        p(1:3) = p(1:3)./data_aspect_ratio;
        p = R_el*R_az*p;
        p = T_cam*p;
        pnt_cam(:,kk) = M_orth*p;
    end
    x_vec = pnt_cam(:,2)-pnt_cam(:,1);
    y_vec = pnt_cam(:,4)-pnt_cam(:,1);
    z_vec = pnt_cam(:,5)-pnt_cam(:,1);
    % find the most-left point
    [~,ix] = min(pnt_cam(1,:));
    % find the lowest point
    [~,iy] = min(pnt_cam(2,:));
    switch ix % use ix to determine z axis, then use iy to determine x and y axis
        case {1,5}
            trans_pnt_z = [x_lim(1),y_lim(1),(z_lim(1)+z_lim(2))/2]';
            if iy == 2
                trans_vec_x_edge = pnt_cam(1:2,2)-pnt_cam(1:2,3);
                trans_vec_y_edge = pnt_cam(1:2,2)-pnt_cam(1:2,1);
                trans_pnt_x = (pnt(1:3,2)+pnt(1:3,1))/2;
                trans_pnt_y = (pnt(1:3,3)+pnt(1:3,2))/2;
            else % iy == 4
                trans_vec_x_edge = pnt_cam(1:2,4)-pnt_cam(1:2,1);
                trans_vec_y_edge = pnt_cam(1:2,4)-pnt_cam(1:2,3);
                trans_pnt_y = (pnt(1:3,4)+pnt(1:3,1))/2;
                trans_pnt_x = (pnt(1:3,3)+pnt(1:3,4))/2;
            end
        case {2,6}
            trans_pnt_z = [x_lim(2),y_lim(1),(z_lim(1)+z_lim(2))/2]';
            if iy == 3
                trans_vec_x_edge = pnt_cam(1:2,3)-pnt_cam(1:2,2);
                trans_vec_y_edge = pnt_cam(1:2,3)-pnt_cam(1:2,4);
                trans_pnt_y = (pnt(1:3,3)+pnt(1:3,2))/2;
                trans_pnt_x = (pnt(1:3,4)+pnt(1:3,3))/2;
            else % iy == 1
                trans_vec_x_edge = pnt_cam(1:2,1)-pnt_cam(1:2,4);
                trans_vec_y_edge = pnt_cam(1:2,1)-pnt_cam(1:2,2);
                trans_pnt_x = (pnt(1:3,1)+pnt(1:3,2))/2;
                trans_pnt_y = (pnt(1:3,4)+pnt(1:3,1))/2;
            end
        case {3,7}
            trans_pnt_z = [x_lim(2),y_lim(2),(z_lim(1)+z_lim(2))/2]';
            if iy == 4
                trans_vec_x_edge = pnt_cam(1:2,4)-pnt_cam(1:2,1);
                trans_vec_y_edge = pnt_cam(1:2,4)-pnt_cam(1:2,3);
                trans_pnt_x = (pnt(1:3,4)+pnt(1:3,3))/2;
                trans_pnt_y = (pnt(1:3,1)+pnt(1:3,4))/2;
            else % iy == 2
                trans_vec_x_edge = pnt_cam(1:2,2)-pnt_cam(1:2,3);
                trans_vec_y_edge = pnt_cam(1:2,2)-pnt_cam(1:2,1);
                trans_pnt_y = (pnt(1:3,2)+pnt(1:3,3))/2;
                trans_pnt_x = (pnt(1:3,1)+pnt(1:3,2))/2;
            end
        case {4,8}
            trans_pnt_z = [x_lim(1),y_lim(2),(z_lim(1)+z_lim(2))/2]';
            if iy == 1
                trans_vec_x_edge = pnt_cam(1:2,1)-pnt_cam(1:2,4);
                trans_vec_y_edge = pnt_cam(1:2,1)-pnt_cam(1:2,2);
                trans_pnt_y = (pnt(1:3,1)+pnt(1:3,4))/2;
                trans_pnt_x = (pnt(1:3,2)+pnt(1:3,1))/2;
            else % iy == 3
                trans_vec_x_edge = pnt_cam(1:2,3)-pnt_cam(1:2,2);
                trans_vec_y_edge = pnt_cam(1:2,3)-pnt_cam(1:2,4);
                trans_pnt_x = (pnt(1:3,3)+pnt(1:3,4))/2;
                trans_pnt_y = (pnt(1:3,2)+pnt(1:3,3))/2;
            end
    end
else
    for kk = 1:8
        p = pnt(:,kk);
        p = T_tar*p;
        p(1:3) = p(1:3)./data_aspect_ratio;
        p = R_el*R_az*p;
        p = T_cam*p;
        z = p(3);
        M_pers = [-N/z,0,0,0;0,-N/z,0,0;0,0,-N,0;0,0,0,1];
        pnt_cam(:,kk) = M_pers*p;
    end
    % find the most-left point
    [~,ix] = min(pnt_cam(1,:));
    % find the lowest point
    [~,iy] = min(pnt_cam(2,:));
    % find the proper axis to label
    switch ix % use ix to determine z axis, then use iy to determine x and y axis
        case {1,5}
            z_vec = pnt_cam(:,5)-pnt_cam(:,1);
            trans_pnt_z = [x_lim(1),y_lim(1),(z_lim(1)+z_lim(2))/2]';
            if iy == 2
                x_vec = pnt_cam(:,2)-pnt_cam(:,1);
                y_vec = pnt_cam(:,3)-pnt_cam(:,2);
                trans_vec_x_edge = pnt_cam(1:2,2)-pnt_cam(1:2,3);
                trans_vec_y_edge = pnt_cam(1:2,2)-pnt_cam(1:2,1);
                trans_pnt_x = (pnt(1:3,2)+pnt(1:3,1))/2;
                trans_pnt_y = (pnt(1:3,3)+pnt(1:3,2))/2;
            else % iy == 4
                y_vec = pnt_cam(:,4)-pnt_cam(:,1);
                x_vec = pnt_cam(:,3)-pnt_cam(:,4);
                trans_vec_x_edge = pnt_cam(1:2,4)-pnt_cam(1:2,1);
                trans_vec_y_edge = pnt_cam(1:2,4)-pnt_cam(1:2,3);
                trans_pnt_y = (pnt(1:3,4)+pnt(1:3,1))/2;
                trans_pnt_x = (pnt(1:3,3)+pnt(1:3,4))/2;
            end
        case {2,6}
            z_vec = pnt_cam(:,6)-pnt_cam(:,2);
            trans_pnt_z = [x_lim(2),y_lim(1),(z_lim(1)+z_lim(2))/2]';
            if iy == 3
                y_vec = pnt_cam(:,3)-pnt_cam(:,2);
                x_vec = pnt_cam(:,4)-pnt_cam(:,3);
                trans_vec_x_edge = pnt_cam(1:2,3)-pnt_cam(1:2,2);
                trans_vec_y_edge = pnt_cam(1:2,3)-pnt_cam(1:2,4);
                trans_pnt_y = (pnt(1:3,3)+pnt(1:3,2))/2;
                trans_pnt_x = (pnt(1:3,4)+pnt(1:3,3))/2;
            else % iy == 1
                x_vec = pnt_cam(:,1)-pnt_cam(:,2);
                y_vec = pnt_cam(:,4)-pnt_cam(:,1);
                trans_vec_x_edge = pnt_cam(1:2,1)-pnt_cam(1:2,4);
                trans_vec_y_edge = pnt_cam(1:2,1)-pnt_cam(1:2,2);
                trans_pnt_x = (pnt(1:3,1)+pnt(1:3,2))/2;
                trans_pnt_y = (pnt(1:3,4)+pnt(1:3,1))/2;
            end
        case {3,7}
            z_vec = pnt_cam(:,7)-pnt_cam(:,3);
            trans_pnt_z = [x_lim(2),y_lim(2),(z_lim(1)+z_lim(2))/2]';
            if iy == 4
                x_vec = pnt_cam(:,4)-pnt_cam(:,3);
                y_vec = pnt_cam(:,1)-pnt_cam(:,4);
                trans_vec_x_edge = pnt_cam(1:2,4)-pnt_cam(1:2,1);
                trans_vec_y_edge = pnt_cam(1:2,4)-pnt_cam(1:2,3);
                trans_pnt_x = (pnt(1:3,4)+pnt(1:3,3))/2;
                trans_pnt_y = (pnt(1:3,1)+pnt(1:3,4))/2;
            else % iy == 2
                y_vec = pnt_cam(:,2)-pnt_cam(:,3);
                x_vec = pnt_cam(:,1)-pnt_cam(:,2);
                trans_vec_x_edge = pnt_cam(1:2,2)-pnt_cam(1:2,3);
                trans_vec_y_edge = pnt_cam(1:2,2)-pnt_cam(1:2,1);
                trans_pnt_y = (pnt(1:3,2)+pnt(1:3,3))/2;
                trans_pnt_x = (pnt(1:3,1)+pnt(1:3,2))/2;
            end
        case {4,8}
            z_vec = pnt_cam(:,8)-pnt_cam(:,4);
            trans_pnt_z = [x_lim(1),y_lim(2),(z_lim(1)+z_lim(2))/2]';
            if iy == 1
                y_vec = pnt_cam(:,1)-pnt_cam(:,4);
                x_vec = pnt_cam(:,2)-pnt_cam(:,1);
                trans_vec_x_edge = pnt_cam(1:2,1)-pnt_cam(1:2,4);
                trans_vec_y_edge = pnt_cam(1:2,1)-pnt_cam(1:2,2);
                trans_pnt_y = (pnt(1:3,1)+pnt(1:3,4))/2;
                trans_pnt_x = (pnt(1:3,2)+pnt(1:3,1))/2;
            else % iy == 3
                x_vec = pnt_cam(:,3)-pnt_cam(:,4);
                y_vec = pnt_cam(:,2)-pnt_cam(:,3);
                trans_vec_x_edge = pnt_cam(1:2,3)-pnt_cam(1:2,2);
                trans_vec_y_edge = pnt_cam(1:2,3)-pnt_cam(1:2,4);
                trans_pnt_x = (pnt(1:3,3)+pnt(1:3,4))/2;
                trans_pnt_y = (pnt(1:3,2)+pnt(1:3,3))/2;
            end
    end
end

if trans_mode == 1
    trans_vec_x = [x_vec(2),-x_vec(1)]';
    trans_vec_y = [y_vec(2),-y_vec(1)]';
    trans_vec_z = [-z_vec(2),z_vec(1)]';
    if trans_vec_x'*trans_vec_x_edge < 0
        trans_vec_x = -trans_vec_x;
    end
    if trans_vec_y'*trans_vec_y_edge < 0
        trans_vec_y = -trans_vec_y;
    end
elseif trans_mode == 2
    trans_vec_x = trans_vec_x_edge;
    trans_vec_y = trans_vec_y_edge;
    trans_vec_z = [-z_vec(2),z_vec(1)]';
end

% Normalize translation vectors
trans_vec_x = trans_vec_x/norm(trans_vec_x,2);
trans_vec_y = trans_vec_y/norm(trans_vec_y,2);
trans_vec_z = trans_vec_z/norm(trans_vec_z,2);

% Compute rotation angles
theta_x = atan2d(x_vec(2),x_vec(1));
theta_y = atan2d(y_vec(2),y_vec(1));
theta_z = atan2d(z_vec(2),z_vec(1));
if abs(theta_x) >= 90
    theta_x = theta_x+180;
end
if abs(theta_y) >= 90
    theta_y = theta_y+180;
end
% if abs(theta_z) >= 90
%     theta_z = theta_z+180;
% end
set(get(h_a,'xlabel'),'rotation',theta_x);
set(get(h_a,'ylabel'),'rotation',theta_y);
set(get(h_a,'zlabel'),'rotation',theta_z);

%% Axis label translation
axislabel_translation(h_a,...
    [trans_pnt_x,trans_pnt_y,trans_pnt_z],...
    [trans_vec_x,trans_vec_y,trans_vec_z],...
    'character')
end
