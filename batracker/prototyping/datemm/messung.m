% Mess-Reihe zu DATEMM
%   She, 26.05.06

global DMSG
DMSG = 0;

Fs = 96000;
vs = 334;
blocklength = 4096;
duration = 10;

path='/scratch/messung/loudspkr/';

% Raum-Referenz (wird nicht benoetigt)
p1 = [0 0 0]'; p2 = [3.52 0 0]'; p3 = [0 2 0]'; p4 = [3.52 2 0]'; p5 = [0 0 2]'; p6 = [3.52 0 2]'; p7 = [0 2 2]'; p8 = [3.52 2 2]';

inside = [1 1 1]';
wall(1).pos1 = p1; wall(1).pos2 = p2; wall(1).pos3 = p4; wall(1).pos4 = p3; 
wall(2).pos1 = p1; wall(2).pos2 = p2; wall(2).pos3 = p6; wall(2).pos4 = p5; 
wall(3).pos1 = p1; wall(3).pos2 = p3; wall(3).pos3 = p7; wall(3).pos4 = p5; 
wall(4).pos1 = p2; wall(4).pos2 = p4; wall(4).pos3 = p8; wall(4).pos4 = p6; 
wall(5).pos1 = p3; wall(5).pos2 = p4; wall(5).pos3 = p8; wall(5).pos4 = p7; 
wall(6).pos1 = p5; wall(6).pos2 = p6; wall(6).pos3 = p8; wall(6).pos4 = p7; 

for w = 1:length(wall)
  dummy = cross(wall(w).pos2-wall(w).pos1,wall(w).pos2-wall(w).pos3);
  wall(w).nv = dummy/norm(dummy)*sign(dot(wall(w).pos1-inside,dummy));
end
clear p1 p2 p3 p4 p5 p6 p7 p8 dummy inside w gamma

src(1).pos=[1.6,.18,.69]';     src(1).name='S1';
src(2).pos=[2.96,1.11,1.09]';  src(2).name='S2';

mic(1).pos=[1.565,1.58,.90]';  mic(1).name = 'M1';
mic(2).pos=[2.10,.32,1.54]';   mic(2).name = 'M2';
mic(3).pos=[.43,.43,.65]';     mic(3).name = 'M3';
mic(4).pos=[2.12,.31,.765]';   mic(4).name = 'M4';
mic(5).pos=[2.715,.22,1.68]';  mic(5).name = 'M5';
mic(6).pos=[2.58,1.74,1.4]';   mic(6).name = 'M6';
mic(7).pos=[.71,1.555,1.655]'; mic(7).name = 'M7';
mic(8).pos=[2.13,1.035,1.61]'; mic(8).name = 'M8';

for block=[0:Fs*duration/blocklength*2]
  start=block*blocklength/2;
  bereich=[start+1,start+blocklength];
  tmp = wavread([path,'a1.wav'],bereich); mic(1).wav = tmp(:,1);
  tmp = wavread([path,'a2.wav'],bereich); mic(2).wav = tmp(:,1);
  tmp = wavread([path,'a3.wav'],bereich); mic(3).wav = tmp(:,1);
  tmp = wavread([path,'a4.wav'],bereich); mic(4).wav = tmp(:,1);
  tmp = wavread([path,'a5.wav'],bereich); mic(5).wav = tmp(:,1);
  tmp = wavread([path,'a6.wav'],bereich); mic(6).wav = tmp(:,1);
  tmp = wavread([path,'a7.wav'],bereich); mic(7).wav = tmp(:,1);
  tmp = wavread([path,'a8.wav'],bereich); mic(8).wav = tmp(:,1);
  main_18
  if length(istqual)>1
    save([path,'datemm',num2str(block)]);
  end

  clear CTRIP DP VTRIP gsrc mix pair trip;
end

% Schaubild mit allen Lokalisierungspunkten

clf; hold on;
cd (path);
logx=dir('*.mat');
for axx=1:length(logx)
  clear istpos;
  load(logx(axx).name);
  for b=1:size(istpos,2)
    if (min(istpos(:,b)>-.5) & (max(istpos(:,b))<4))
      plot3(istpos(1,b),istpos(2,b),istpos(3,b),'.');
    end
  end
end
