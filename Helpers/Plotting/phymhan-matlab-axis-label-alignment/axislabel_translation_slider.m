function axislabel_translation_slider
% Create UI figure window and components

global AXISALIGN_TRANS_A AXISALIGN_TRANS_B
if isempty(AXISALIGN_TRANS_A), AXISALIGN_TRANS_A = 1; end
if isempty(AXISALIGN_TRANS_B), AXISALIGN_TRANS_B = 1; end

h_f = uifigure('Position',[100 100 350 275]);

% cg = uigauge(fig,'Position',[100 100 120 120]);

sld1 = uislider(h_f,...
    'Position',[50 200 250 3],...
    'Limits',[0 5],...
    'ValueChangingFcn',@(sld,event) sliderMoving1(event));
sld1.Value = 1;
sld2 = uislider(h_f,...
    'Position',[50 75 250 3],...
    'Limits',[0 5],...
    'ValueChangingFcn',@(sld,event) sliderMoving2(event));
sld2.Value = 1;
end

% Create ValueChangingFcn callback
function sliderMoving1(event)
global AXISALIGN_TRANS_A
AXISALIGN_TRANS_A = event.Value;
align_axislabel([],gca);
end

function sliderMoving2(event)
global AXISALIGN_TRANS_B
AXISALIGN_TRANS_B = event.Value;
align_axislabel([],gca);
end
