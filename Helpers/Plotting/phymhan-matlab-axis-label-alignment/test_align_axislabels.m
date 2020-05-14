function test_align_axislabels(~,~)
h_a = gca;
[thx,thy,thz] = axislabel_rotation_angle(h_a);
set(get(h_a,'xlabel'),'rotation',thx(1));
set(get(h_a,'ylabel'),'rotation',thy(1));
set(get(h_a,'zlabel'),'rotation',thz(1));
end
