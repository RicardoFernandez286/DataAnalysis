% pnt = [-2,-2,-4,1;-2,2,-4,1;2,2,-4,1;2,-2,-4,1;...
%     -2,-2,-8,1;-2,2,-8,1;2,2,-8,1;2,-2,-8,1]';
% F = 200;
% N = 40;
% a = -(F+N)/(F-N);
% b = -2*F*N/(F-N);
% thx = 0;
% thy = 0;
% thz = 0;
% 
% for k = 1:1
% thx = 1/180*pi;
% thy = thy+1/180*pi;
% Rx = [1,0,0,0;0,cos(-thx),-sin(-thx),0;0,sin(-thx),cos(-thx),0;0,0,0,1];
% Ry = [cos(-thy),0,sin(-thy),0;0,1,0,0;-sin(-thy),0,cos(-thy),0;0,0,0,1];
% Rz = [cos(-thz),-sin(-thz),0,0;sin(-thz),cos(-thz),0,0;0,0,1,0;0,0,0,1];
% pnt_r = Rx*Ry*Rz*pnt;
% %pnt_cam = M_verl*pnt_r;
% pnt_cam = zeros(size(pnt));
% T = [1,0,0,0;0,1,0,0;0,0,1,N-6;0,0,0,1];
% for kk = 1:size(pnt,2)
%     p = pnt_r(:,kk);
%     %pnt_cam(:,kk) = M_in*T*p;
%     p = T*p;
%     z = p(3);
%     M_in = [-N/z,0,0,0;0,-N/z,0,0;0,0,-N,0;0,0,0,1];
%     pnt_cam(:,kk) = M_in*p;
%     
% end
% pnt_cam
% cla
% hold on
% plot(pnt_cam(1,[1:4,1]),pnt_cam(2,[1:4,1]),'r')
% plot(pnt_cam(1,[5:8,5]),pnt_cam(2,[5:8,5]),'b')
% plot(pnt_cam(1,[1,5]),pnt_cam(2,[1,5]),'g')
% plot(pnt_cam(1,[2,6]),pnt_cam(2,[2,6]),'g')
% plot(pnt_cam(1,[3,7]),pnt_cam(2,[3,7]),'g')
% plot(pnt_cam(1,[4,8]),pnt_cam(2,[4,8]),'g')
% pause(0.1)
% end
