function [theta_x,theta_y,theta_z] = axislabel_rotation_angle(h_a)
% This function computes the proper angles of x, y and z labels in order to
% align them to the direction of their axes.
% It can be used as part of an 'ActionPostCallback', which is a callback 
% function of a rotate3d object.
% Each of the output arguments theta_x, theta_y and theta_z is a 1-by-4 
% vector whose elements represent the angles of the four edges of the axes 
% box. In each group, the four edges are listed anticlockwise.
% The function still works when the projection mode is perspective or when 
% the data aspect ratio is not [1 1 1].
% Compared with 'axislabel_rotation', the process that determines which edge 
% is the proper tick-marked axis is removed and all of the twelve angles 
% are returned.
%
% Feb/05/2015, hanligong@gmail.com

proj_mode = get(h_a,'projection');
if strcmpi(proj_mode,'orthographic')
    proj_mode = 1;
else
    proj_mode = 2;
end
data_aspect_ratio = get(h_a,'dataaspectratio');
data_aspect_ratio = data_aspect_ratio(:);
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
    for kk = 1:8
        p = pnt(:,kk);
        p = T_tar*p;
        p(1:3) = p(1:3)./data_aspect_ratio;
        p = R_el*R_az*p;
        p = T_cam*p;
        pnt_cam(:,kk) = M_orth*p;
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
end
x_vec = [pnt_cam(:,2)-pnt_cam(:,1),...
    pnt_cam(:,3)-pnt_cam(:,4),...
    pnt_cam(:,7)-pnt_cam(:,8),...
    pnt_cam(:,6)-pnt_cam(:,5)];
y_vec = [pnt_cam(:,4)-pnt_cam(:,1),...
    pnt_cam(:,3)-pnt_cam(:,2),...
    pnt_cam(:,7)-pnt_cam(:,6),...
    pnt_cam(:,8)-pnt_cam(:,5)];
z_vec = [pnt_cam(:,5)-pnt_cam(:,1),...
    pnt_cam(:,6)-pnt_cam(:,2),...
    pnt_cam(:,7)-pnt_cam(:,3),...
    pnt_cam(:,8)-pnt_cam(:,4)];
theta_x = atan2d(x_vec(2,:),x_vec(1,:));
theta_y = atan2d(y_vec(2,:),y_vec(1,:));
theta_z = atan2d(z_vec(2,:),z_vec(1,:));
end
