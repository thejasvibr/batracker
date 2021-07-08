% FUNCTION plotwall(wall[,fig])
%  - akustisches Raummodell - She, 26.07.05
%
% zeichnet Waende wall in 3D-Raum fig ein

function plotwall(wall,fig)

if nargin < 2 fig = 0; end
if fig ~= 0 figure(fig); end
axis equal; hold on; box on;

for w = 1:length(wall)
% Seitenflaeche
%   edges=[wall(w).pos1,wall(w).pos2,wall(w).pos3,wall(w).pos4];
%   set(fill3(edges(1,:),edges(2,:),edges(3,:),'k'),'FaceAlpha',.1,'EdgeAlpha',.1);
  
% Kanten
%  l = [wall(w).pos1,wall(w).pos2,wall(w).pos3,wall(w).pos1,wall(w).pos4,wall(w).pos3,wall(w).pos2,wall(w).pos4]';
  l = [wall(w).pos1,wall(w).pos2,wall(w).pos3,wall(w).pos4,wall(w).pos1]';
  line(l(:,1),l(:,2),l(:,3));
  
% Normalenvektor  
  l = [.25*(wall(w).pos1+wall(w).pos2+wall(w).pos3+wall(w).pos4),.25*(wall(w).pos1+wall(w).pos2+wall(w).pos3+wall(w).pos4)+.1*wall(w).nv]';
  plot3(l(:,1),l(:,2),l(:,3),'g');
end

return;
