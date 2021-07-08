% FUNCTION plotsrc(src[,fig])
%  - akustisches Raummodell - She, 03.09.03
%
% zeichnet Position der Schallquellen src im 3D-Raum fig

function plotmic(src,fig)

if nargin < 2 fig = 0; end
if fig ~= 0 figure(fig); end
axis equal; hold on; box on;

for s = 1:length(src)
  if isfield (src(s),'dir')
    d = [src(s).pos,src(s).pos+.1/sqrt(src(s).dir'*src(s).dir)*src(s).dir]';
    line(d(:,1),d(:,2),d(:,3))
  end
  plot3(src(s).pos(1),src(s).pos(2),src(s).pos(3),'*')
end

return;
