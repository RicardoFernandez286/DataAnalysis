function newmap = darkb2r(cmin_input,cmax_input,n_total,n_whites)
%BLUEWHITERED   dark Blue, white, and dark red color map.
%   this matlab file is designed to draw anomaly figures. the color of
%   the colorbar is from dark blue to white and then to dark red, corresponding to 
%   the anomaly values from negative to zero to positive, respectively. 
%   The color white always correspondes to value zero. 
%   
%   You should input two values like caxis in matlab, that is the min and
%   the max value of color values designed.  e.g. colormap(darkb2r(-3,5))
%   
%   the brightness of blue and red will change according to your setting,
%   so that the brightness of the color corresponded to the color of his
%   opposite number
%   e.g. colormap(darkb2r(-3,6))   is from light blue to deep red
%   e.g. colormap(darkb2r(-3,3))   is from deep blue to deep red
%
%   I'd advise you to use colorbar first to make sure the caxis' cmax and cmin.
%   Besides, there is also another similar colorbar named 'b2r', in which the 
%   color is lighter.
%
%   by Cunjie Zhang, 2011-3-14
%   find bugs ====> email : daisy19880411@126.com
%   updated:  Robert Beckman help to fix the bug when start point is zero, 2015-04-08
%   Updated:  Ricardo Fernandez-Teran added "n_whites" to the input options (2018.05.29)
%   Updated:  Ricardo Fernandez-Teran made the function compatible with UIAxes by removing the GCA command (2019.04.07)
%  
%   Examples:
%   ------------------------------
%   figure
%   peaks;
%   colormap(darkb2r(-6,8)), colorbar, title('darkb2r')
%   


%% check the input
if nargin > 4
   warning('Input two variables, the range of caxis , for example : colormap(b2r(-3,3)) or three, b2r(-3,3,ncolors)');
elseif nargin > 2 && nargin < 4
   n_total = 250;
   n_whites = 1;
end

if cmin_input >= cmax_input
    disp('input error');
    disp('the color range must be from a smaller one to a larger one');
end

%% control the figure caxis 
lims = [0 1]; % Get figure caxis formation

%% color configuration : from blue to to white then to red
red_top     = [0.5 0 0];
red_middle  = [1 0 0];
white_middle= [1 1 1];
blue_middle = [0 0 1];
blue_bottom = [0 0 0.5];

if mod(n_whites,2) ~= 0
    n_whites    = n_whites - 1;
end

n_reds      = round((n_total - n_whites)/2)+1;
n_blues     = round((n_total - n_whites)/2)+1;

%% Color interpolation 
m=2/3; % Position of the set colors
color_input     = [blue_bottom; blue_middle; white_middle];
oldsteps        = [-1 -m 0];
newsteps        = linspace(-1, 0, n_blues);  
bluemap         = zeros(n_blues,3);
for j=1:3
   bluemap(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
end

color_input     = [white_middle; red_middle; red_top];
oldsteps        = [0 m 1];
newsteps        = linspace(0, 1, n_reds);  
redmap          = zeros(n_reds,3);
for j=1:3
   redmap(:,j)  = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
end

bluemap(end,:)  = [];
redmap(1,:)     = [];

if n_whites ~= 0
    whitemap        = ones(n_whites,3);
    newmap_all      = [bluemap;whitemap;redmap];
else
    newmap_all      = [bluemap;redmap];
end

%% Category Discussion according to the cmin and cmax input

%  the color data will be remaped to color range from -max(abs(cmin_input),cmax_input)
%  to max(abs(cmin_input),cmax_input) , and then squeeze the color data
%  in order to make sure the blue and red color selected corresponded
%  to their math values

%  for example :
%  if darkb2r(-3,6) ,the color range is from light blue to deep red , so that
%  the light blue valued at -3 correspondes to light red valued at 3


%% Category Discussion according to the cmin and cmax input
% first : from negative to positive
% then  : from positive to positive
% last  : from negative to negative


if (cmin_input < 0)  &&  (cmax_input > 0)  
        
    if abs(cmin_input) < cmax_input 
         
        % |--------|---------|--------------------|    
      % -cmax      cmin       0                  cmax         [cmin,cmax]
 
       start_point = max(round((cmin_input+cmax_input)/2/cmax_input*n_total),1);
       newmap = squeeze(newmap_all(start_point:n_total,:));
       
    elseif abs(cmin_input) >= cmax_input
        
         % |------------------|------|--------------|    
       %  cmin                0     cmax          -cmin         [cmin,cmax]   
       
       end_point = max(round((cmax_input-cmin_input)/2/abs(cmin_input)*n_total),1);
       newmap = squeeze(newmap_all(1:end_point,:));
    end
       
elseif cmin_input >= 0

       if lims(1) < 0 
           disp('caution:')
           disp('there are still values smaller than 0, but cmin is larger than 0.')
           disp('some area will be in red color while it should be in blue color')
       end
        % |-----------------|-------|-------------|    
      % -cmax               0      cmin          cmax         [cmin,cmax]
 
       start_point = max(round((cmin_input+cmax_input)/2/cmax_input*n_total),1);
       newmap = squeeze(newmap_all(start_point:n_total,:));

elseif cmax_input <= 0
       if lims(2) > 0 
           disp('caution:')
           disp('there are still values larger than 0, but cmax is smaller than 0.')
           disp('some area will be in blue color while it should be in red color')
       end
       
         % |------------|------|--------------------|    
       %  cmin         cmax    0                  -cmin         [cmin,cmax]      
       end_point = max(round((cmax_input-cmin_input)/2/abs(cmin_input)*n_total),1);
       newmap = squeeze(newmap_all(1:end_point,:));
end
