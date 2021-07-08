% Raum-Simulation fuer TDOA-Verfahren
%
%   She, 26.07.05

global vs Fs sint DMSG;

DMSG = 1;
INTERACTIVE = 1;

clear mic src wall
mic(1) = struct('pos',[],'name','','wav',[]);
src(1) = struct('pos',[],'rd',0,'name','','wav',[]);
wall(1) = struct('pos1',[],'pos2',[],'pos3',[],'pos4',[],'gamma',1,'name','','nv',[]);

rand('state',sum(100*clock))

if DMSG>=3 echo on; end;

% Parameter
vs = 334;
Fs = 96000;
numsrc = 2;
nummic = 6;
spmorder = 3;
numsamples = 8000;
room = 'room_office';

if INTERACTIVE
  disp('manipulate values, then type return');
  disp(['  numsrc = ', num2str(numsrc), ...
    ', nummic = ', num2str(nummic), ...
    ', Fs = ', num2str(Fs), ...
    ', spmorder = ', num2str(spmorder), ...
    ', numsamples = ',num2str(numsamples), ...
    ', room = ',room]);
keyboard
end

% Waende
if DMSG disp(['... loading walls from ',room]); end;
cd spm
run (room);

% zufaellige Geometrie im Bereich 10%-90% zwischen den Waenden
% mit gegenseitigem Abstand > 1m
if DMSG disp('... locating at random positions inside the room'); end;
wmin = min([wall.pos1,wall.pos2,wall.pos3,wall.pos4]')';
wmax = max([wall.pos1,wall.pos2,wall.pos3,wall.pos4]')';
abstand_ok = 0;
while ~abstand_ok
  abstand_ok = 1;
  pos = repmat(.8*(wmax-wmin),1,numsrc+nummic).*rand(3,numsrc+nummic)+repmat(wmin+.1*wmax,1,numsrc+nummic);
  for n=1:numsrc+nummic-1
    for m=n+1:numsrc+nummic
      if norm(pos(:,n)-pos(:,m))<1
        abstand_ok = 0;
      end
    end
  end
end

for n=1:numsrc
  src(n).pos = pos(:,n);
  src(n).rd = 0;
  src(n).name = ['S',num2str(n)];
end
for n=1:nummic
  mic(n).pos = pos(:,n+numsrc);
  mic(n).name = ['M',num2str(n)];
end

spkdelta=zeros(1,1,numsrc);
clear n m abstand_ok room pos wmin wmax

if INTERACTIVE 
  disp('change mic(:).pos or src(:).pos');
  disp('add src(:).fn for input wav, otherwise white random');
  disp('to view room run plotroom');
  cd ..
  keyboard; 
  cd spm
end;

if DMSG disp('... calculating true delta'); end;
for s=1:length(src)
  for m1=1:length(mic)
    for m2=1:length(mic)
      spkdelta(m1,m2,s) = round((norm(mic(m1).pos-src(s).pos)-norm(mic(m2).pos-src(s).pos))/vs*Fs);
    end
  end
end
maxmicdist = 0;
for m1=1:length(mic)-1
  for m2=m1:length(mic)
    maxmicdist = max([maxmicdist,ceil(norm(mic(m1).pos-mic(m2).pos)*Fs/vs)]);
  end
end

% Spiegelquellen
if DMSG disp('... calculating virtual sources'); end;
srcn = virtsrc(src,wall,spmorder);
if DMSG disp(['total: ',num2str(length(srcn)),' sources']); end;
  
maxsrcdist = 0;
for s=1:length(src)
  for m1=1:length(mic)
    maxsrcdist = max([maxsrcdist,ceil(norm(mic(m1).pos-src(s).pos)*Fs/vs)]);
  end
end

% Sample-Ausschnitt
sint(1) = spmorder*(maxmicdist+maxsrcdist);
sint(2) = sint(1)+numsamples-1;
if DMSG disp(['... interval is [',num2str(sint(1)),',',num2str(sint(2)),']']); end;

if DMSG disp('... setting src(:).wav'); end;
for n=1:numsrc
  src(n).wav = -1+2*rand(sint(2),1);
  if isfield(src(n),'fn')
    if ~isempty(src(n).fn)
      src(n).wav = wavread(src(n).fn,sint(2));
    end
  end
end

% Spiegelquellen
if DMSG disp('... calculating virtual sources'); end;
src = virtsrc(src,wall,spmorder);
if DMSG disp(['total: ',num2str(length(src)),' sources']); end;

% Schallausbreitung
if DMSG disp('... calculating microphone signal'); end;
mic = propag_tr(mic,src,wall);

for n=1:length(src)
  src(n).wav=[];
end

echo off
cd ..
clear m1 m2 maxmicdist n nummic numsamples numsrc s sint spmorder srcn INTERACTIVE DMSG
