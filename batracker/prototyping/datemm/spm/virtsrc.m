% FUNCTION srcn = virtsrc(src,wall,rdmax[,dm])
%  - akustisches Raummodell - She, 26.07.05
%
% Rekusive Berechnung virtueller Akustikquellen mit der Spiegelungs-
% methode aus tatsaechlichen bzw. schon berechneten virtuellen Quellen 
% src und den Waenden wall. Beruecksichtigt werden Spiegelungen an 
% allen Waenden wall ausser denen mit Index dm, wobei die Spiegelung
% hoechstens von der Ordnung rdmax ist.
% Ungueltige Spiegelquellen (entstanden durch Spiegelung an nicht
% reflektierender Wandseite) werden vermieden (Wandnormalenvektor nach aussen).
%
% Bem: Ausgabe enthaelt sowohl echte als auch virtuelle Quellen.
%      Beruecksichtigt Wand-Reflexionskoeffizient falls vorhanden.
%      Default: dm=0

function srcn = virtsrc(src,wall,rdmax,dm)
if nargin < 4 dm = 0; end
  
global dmsg;  

srcn = src;

for s = 1:length(src)
% Rekursionsabbruch falls maximale Tiefe erreicht
  if src(s).rd < rdmax
    for w = 1:length(wall)
% Gueltigkeitspruefung
      if sign(dot(src(s).pos-wall(w).pos1,wall(w).nv)) == -1
        if dmsg>=2 disp (sprintf('spiegle Quelle %s an Wand %s',src(s).name,wall(w).name)); end;
% Variablenstruktur und gleiche Inhalte uebernehmen
        newsrc = src(s);
% Spiegelung src(s) an Ebene wall(w)
        newsrc.pos = src(s).pos+2*dot((wall(w).pos2-src(s).pos),wall(w).nv)*wall(w).nv;
        newsrc.rd = src(s).rd+1;
        newsrc.name = [src(s).name '_' wall(w).name ];
% Reflexionskoeffizient
        if isfield (wall(w),'gamma') newsrc.wav = src(s).wav*wall(w).gamma; end
% Ausrichtung
        if isfield (src(s),'dir')
          newsrc.dir = src(s).dir-2*dot(src(s).dir,wall(w).nv)/dot(wall(w).nv,wall(w).nv)*wall(w).nv;
        end
% Rekursiver Aufruf
        srcn = [srcn virtsrc(newsrc,wall,rdmax,w)];
      end
    end
  end
end

return;
