% FUNCTION plotmic(mic[,fig])
%  - akustisches Raummodell - She, 03.09.03
%
% zeichnet Position und ggf. Ausrichtung der Mikrofone 
% mic im 3D-Raum fig

function plotmic(mic,fig)

if nargin < 2 fig = 0; end
if fig ~= 0 figure(fig); end
axis equal; hold on; box on;

for m = 1:length(mic)
  if isfield (mic(m),'dir')
    d = [mic(m).pos,mic(m).pos+.1/sqrt(mic(m).dir'*mic(m).dir)*mic(m).dir]';
    line(d(:,1),d(:,2),d(:,3))
  end
  plot3(mic(m).pos(1),mic(m).pos(2),mic(m).pos(3),'o')
end

return;
