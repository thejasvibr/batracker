% echo chamber at LSS
%   She, 26.07.05

% edge-points
p1 = [0 0 0]';
p2 = [4.99 0 0]';
p3 = [0 5.63 0]';
p4 = [4.11 4.86 0]';
p5 = [0 0 5.54]';
p6 = [4.99 0 5.00]';
p7 = [0 5.63 5.54]';
p8 = [4.11 4.86 5.14]';

% reflection coefficient
gamma = 0.9;

% point inside
inside = [1 1 1]';

wall(1).pos1 = p1; wall(1).pos2 = p2; wall(1).pos3 = p4; wall(1).pos4 = p3; wall(1).gamma = gamma; wall(1).name = 'W1';
wall(2).pos1 = p1; wall(2).pos2 = p2; wall(2).pos3 = p6; wall(2).pos4 = p5; wall(2).gamma = gamma; wall(2).name = 'W2';
wall(3).pos1 = p1; wall(3).pos2 = p3; wall(3).pos3 = p7; wall(3).pos4 = p5; wall(3).gamma = gamma; wall(3).name = 'W3';
wall(4).pos1 = p2; wall(4).pos2 = p4; wall(4).pos3 = p8; wall(4).pos4 = p6; wall(4).gamma = gamma; wall(4).name = 'W4';
wall(5).pos1 = p3; wall(5).pos2 = p4; wall(5).pos3 = p8; wall(5).pos4 = p7; wall(5).gamma = gamma; wall(5).name = 'W5';
wall(6).pos1 = p5; wall(6).pos2 = p6; wall(6).pos3 = p8; wall(6).pos4 = p7; wall(6).gamma = gamma; wall(6).name = 'W6';

for w = 1:length(wall)
  dummy = cross(wall(w).pos2-wall(w).pos1,wall(w).pos2-wall(w).pos3);
  wall(w).nv = dummy/norm(dummy)*sign(dot(wall(w).pos1-inside,dummy));
end

clear p1 p2 p3 p4 p5 p6 p7 p8 dummy inside w gamma
