% FUNCTION plotlabel(obj[,fig])
%  - akustisches Raummodell - She, 03.09.03
%
% bezeichnet Objekte obj im 3D-Plot
% obj in { src, mic, wall }

function plotlabel(obj,fig)

if nargin < 2 fig = 0; end
if fig ~= 0 figure(fig); end
hold on;

for o = 1:length(obj)
  text(obj(o).pos(1),obj(o).pos(2),obj(o).pos(3),obj(o).name)
end

return;
