% M-FILE plotroom
%  - akustisches Raummodell - She, 26.07.05
%
% zeichnet Raum und alle Objekte

clf; 
axis equal; hold on; box on;

for m = 1:length(mic)
  if isfield (mic(m),'dir')
    d = [mic(m).pos,mic(m).pos+.1/sqrt(mic(m).dir'*mic(m).dir)*mic(m).dir]';
    line(d(:,1),d(:,2),d(:,3))
  end
  plot3(mic(m).pos(1),mic(m).pos(2),mic(m).pos(3),'o')
  text(mic(m).pos(1),mic(m).pos(2),mic(m).pos(3),mic(m).name)
%    d = [mic(m).pos,[mic(m).pos(1),mic(m).pos(2),0]']';
%    line(d(:,1),d(:,2),d(:,3))
end
for s = 1:length(src)
  if isfield (src(s),'dir')
    d = [src(s).pos,src(s).pos+.1/sqrt(src(s).dir'*src(s).dir)*src(s).dir]';
    line(d(:,1),d(:,2),d(:,3))
  end
  plot3(src(s).pos(1),src(s).pos(2),src(s).pos(3),'*')
  text(src(s).pos(1),src(s).pos(2),src(s).pos(3),src(s).name)
%    d = [src(s).pos,[src(s).pos(1),src(s).pos(2),0]']';
%    line(d(:,1),d(:,2),d(:,3))
end
for w = 1:length(wall)
  l = [wall(w).pos1,wall(w).pos2,wall(w).pos3,wall(w).pos4,wall(w).pos1]';
  line(l(:,1),l(:,2),l(:,3));
%    l = [.25*(wall(w).pos1+wall(w).pos2+wall(w).pos3+wall(w).pos4),.25*(wall(w).pos1+wall(w).pos2+wall(w).pos3+wall(w).pos4)+.1*wall(w).nv]';
%    plot3(l(:,1),l(:,2),l(:,3),'g');
end

clear m s w l d
