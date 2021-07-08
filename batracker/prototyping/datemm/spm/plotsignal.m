% FUNCTION plotsignal(obj[,fig])
%  - akustisches Raummodell - She, 03.09.03
%
% stellt Zeitsignal der Objekte obj dar
% obj in {src, mic }

function plotsignal(obj,fig)

global Fs;

if nargin < 2 fig = 0; end
if fig ~= 0 figure(fig); end
clf
hold on

color = 'bgrcmyk';

for o = 1:length(obj)
  t = 1/Fs*[0:length(obj(o).wav)-1];
  plot(t,obj(o).wav,color(mod(o,7)+1));
end

legend(obj.name);

return;
