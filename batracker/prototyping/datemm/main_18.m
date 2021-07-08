% DATEMM-Verfahren
%   She, 26.05.06

global DMSG;

if (length(who('DMSG'))==0) DMSG = 1; end
if (length(who('INTERACTIVE'))==0) INTERACTIVE = 0; end

if INTERACTIVE echo on; end;

% Stellschrauben
DP.FN = 'dp12';
DP.P.PF = [1 .8 .5 .3];             % Glaettungsfunktion (symmetrisch)
DP.P.AUF = 1;                       % Aufwertgewichtung
DP.P.AB = -1.5;                     % Abwertgewichtung
DP.P.FIND = 'findmaxi';
DP.NA = 10;                         % Maximalzahl Maxima AC
DP.P.NX = 10;                       % Maximalzahl Maxima XC
DP.TUA = 0.3;                       % upper threshold AC
DP.P.TUX = 0.5;                     % upper threshold XC
DP.TLA = .1;                        % lower threshold AC
DP.P.TLX = .3;                      % lower threshold XC
DP.P.DIST = 2;                      % neighbour maximum combination distance
VTRIP.FN = 'vtrip5';
VTRIP.P.PF = [.3 .6 .8 1 .8 .6 .3]; % Glaettungsfunktion (symmetrisch)
CTRIP.FN = 'ctrip4';
ACPHATTH = 4;
PHATTH = 1;
TRIPREL = .2;

% manipulate function calls, then type return
if INTERACTIVE echo off; keyboard; end


clear pair trip gsrc
pair(1) = struct('m1',0,'m2',0,'corr',[]);
trip(1) = struct('mic',0,'pair',0,'index',0,'tdoa',0,'qual',0,'used',0);
gsrc(1) = struct('index',0,'tdoa',0,'tdoaq',0,'valid',0,'qual',0);

N = length(mic(1).wav);
M = length(mic);
MX = floor(.5*M*(M-1));

if DMSG disp ('Auto- und Kreuzkorrelation'); end;
acs = zeros(2*N,1);
for m1=1:M
  mic(m1).F = fft([mic(m1).wav;zeros(N,1)]);
  mic(m1).F2 = mic(m1).F.*conj(mic(m1).F);
  acs = acs+mic(m1).F2;
end
acs=min([ACPHATTH*M*ones(size(mic(1).F))';1./acs'])';
for m1=1:M
  dummy = real(ifft(mic(m1).F2.*acs));
  mic(m1).ac = dummy(1:end/2)./[N:-1:1]'.*[ones(N/2,1);zeros(N/2,1)];
end

mx = 1;
for m1=1:M-1
  for m2=m1+1:M
    dummy1 = mic(m1).F.*conj(mic(m2).F);
    dummy = abs(ifft(dummy1.*min([PHATTH*ones(size(dummy1))';1./abs(dummy1')])'));
    pair(mx).corr = [dummy(end/2+2:end);dummy(1:end/2)]./[[1:N-1],[N:-1:1]]'.*[zeros(N/2,1);ones(N-1,1);zeros(N/2,1)];
    mx = mx+1;
  end
end

clear mx m1 m2 dummy dummy1
for m1=1:M
  mic(m1).F = [];
end


if DMSG disp ('Direktpfad-Erkennung'); end;
% relevante Betragsmaxima der Autokorrelierten  
for m1=1:M
  mic(m1).acm = double(feval(DP.P.FIND,abs(mic(m1).ac(2:end)),DP.NA,DP.TUA,DP.TLA))+1; 
  mic(m1).acmv = abs(mic(m1).ac(mic(m1).acm+1));
end
pi = 0;
mics2pair = zeros(M,M);
for m1=1:M-1
  for m2=m1+1:M
    pi = pi+1;
    pair(pi).m1 = m1;
    pair(pi).m2 = m2;
    mdist = norm(mic(m1).pos-mic(m2).pos)*Fs/vs;

    [pair(pi).rcm,pair(pi).rcmv] = feval(DP.FN,pair(pi).corr,mic(m1).acm,mic(m1).acmv,mic(m2).acm,mic(m2).acmv,mdist,DP.P);
    if DMSG disp (['  Paar ',num2str(m1),'-',num2str(m2),': ',num2str(length(pair(pi).rcm)),' Max']); end;
    if DMSG>=2 disp (['    dp: ',num2str([pair(pi).rcm]')]); disp(' '); end;
    mics2pair(m1,m2)=pi;
    mics2pair(m2,m1)=-pi;
  end
end
clear m1 m2 mdist pi


if DMSG disp ('Tripel-Erkennung'); end;
ti = 0;
for m1=1:M-2
  for m2=m1+1:M-1
    for m3=m2+1:M
      p12 = mics2pair(m1,m2);
      p23 = mics2pair(m2,m3);
      p13 = mics2pair(m1,m3);

      dummy = feval(VTRIP.FN,pair(p12).rcm,pair(p23).rcm,pair(p13).rcm,VTRIP.P);
      if DMSG disp (['  Tripel ',num2str(m1),'-',num2str(m2),'-',num2str(m3),': ',num2str(size(dummy,1)),' Treffer']); end;

      for i=1:size(dummy,1) % alle erkannten TDOA-Tripel dieses Sensor-Tripels
        ti = ti+1;
        trip(ti).mic = [m1,m2,m3]';
        trip(ti).pair = [p12,p23,p13]';
        trip(ti).index = dummy(i,[1:3])';
        trip(ti).tdoa = [pair(p12).rcm(dummy(i,1)),pair(p23).rcm(dummy(i,2)),pair(p13).rcm(dummy(i,3))]';
        trip(ti).qual = (pair(p12).rcmv(dummy(i,1))+pair(p23).rcmv(dummy(i,2))+pair(p13).rcmv(dummy(i,3)))*dummy(i,4);
        trip(ti).used = 0;
      end
    end
  end
end
clear m1 m2 m3 dummy p12 p23 p13 ti i


if DMSG disp ('Tripel zu Zielen kombinieren'); end;
% Sortierung nach absteigendem Qualitaetsmass
[h1,h2] = sortrows([trip(:).qual]');
tr = flipud(h2);

if length(tr)<3
  disp ('Zu wenige Tripel gefunden');
  return;
end

% relevante Tripel
h3 = find(flipud(h1)>TRIPREL*h1(end));
gi = 0;
for i=1:h3(end)
  if ~trip(tr(i)).used
% Zielinitialisierung mit staerkstem Tripel    
    gi = gi+1;
    gsrc(gi).index = zeros(MX,1);
    gsrc(gi).index(trip(tr(i)).pair) = trip(tr(i)).index;
    trip(tr(i)).used = i;
    gsrc(gi).qual = trip(tr(i)).qual;

% Suche nach passenden Kombinationen
    comb1 = double(feval(CTRIP.FN,gsrc(gi).index,[trip(tr(i+1:end)).pair],[trip(tr(i+1:end)).index]));
    if length(comb1)>0
      gsrc(gi).index([trip(tr(comb1+i)).pair]) = [trip(tr(comb1+i)).index];
      gsrc(gi).qual = sum([trip(tr(comb1+i)).qual]);

% zweiter Suchlauf mit erweitertem Datensatz
      comb2 = double(feval(CTRIP.FN,gsrc(gi).index,[trip(tr(i+1:end)).pair],[trip(tr(i+1:end)).index]));
      for h4=1:length(comb1)
        trip(tr(comb1(h4)+i)).used = i;
      end
      if length(comb2)>0
        gsrc(gi).index([trip(tr(comb2+i)).pair]) = [trip(tr(comb2+i)).index];
        gsrc(gi).qual = gsrc(gi).qual+sum([trip(tr(comb2+i)).qual]);
        for h4=1:length(comb2)
          trip(tr(comb2(h4)+i)).used = i;
        end
      end
    end

% Zuordnung Indizes -> Laufzeitfdifferenzen
    pi = find(gsrc(gi).index~=0);
    gsrc(gi).valid = zeros(MX,1);
    gsrc(gi).valid(pi) = 1;
    gsrc(gi).tdoa = zeros(MX,1);
    gsrc(gi).tdoaq = zeros(MX,1);
    for h4=1:length(pi)
      gsrc(gi).tdoa(pi(h4)) = pair(pi(h4)).rcm(gsrc(gi).index(pi(h4)));
      gsrc(gi).tdoaq(pi(h4)) = pair(pi(h4)).rcmv(gsrc(gi).index(pi(h4)));
    end
  end
end
clear h1 h2 h3 h4 i comb1 comb2 pi gi tr nmax

% kombiniere Datensätze

n1 = 1;
while n1<=length(gsrc)-1
  n2 = n1+1;
  while n2<=length(gsrc)
    e = sum(gsrc(n1).valid.*gsrc(n2).valid);
    if e>0 & sum(abs(gsrc(n1).valid.*gsrc(n2).valid.*(gsrc(n1).tdoa-gsrc(n2).tdoa)))<2*e
      c = find(gsrc(n1).valid & gsrc(n2).valid);
      n = find(gsrc(n2).valid-(gsrc(n1).valid & gsrc(n2).valid));
      gsrc(n1).tdoa(c) = (gsrc(n1).tdoa(c)*gsrc(n1).qual+gsrc(n2).tdoa(c)*gsrc(n2).qual)/(gsrc(n1).qual+gsrc(n2).qual);
      gsrc(n1).tdoaq(c) = (gsrc(n1).tdoaq(c)*gsrc(n1).qual+gsrc(n2).tdoaq(c)*gsrc(n2).qual)/(gsrc(n1).qual+gsrc(n2).qual);
      gsrc(n1).tdoa(n) = gsrc(n2).tdoa(n);
      gsrc(n1).tdoaq(n) = gsrc(n2).tdoaq(n);
      gsrc(n1).valid = gsrc(n1).valid | gsrc(n2).valid;
      gsrc(n1).qual = gsrc(n1).qual+gsrc(n2).qual;
      gsrc(n2)=[];
    else
      n2 = n2+1;
    end
  end
  n1 = n1+1;
end

[h1,h2] = sortrows([gsrc.qual]');
ord = flipud(h2);

clear h1 h2 n1 n2 e c n

echo off

% Kontrollausgabe
if (length(who('soll'))>0)
  soll = [];
  for i=1:M-1
    soll = [soll;shiftdim(spkdelta(i,[i+1:M],:),1)];
  end
  soll

  ist = [gsrc(ord).tdoa;gsrc(ord).qual]

  sollpos = [src(1:size(spkdelta,3)).pos]
end

istpos = [];
istqual = [];
n = 1;
for i=1:length(gsrc)
  if sum(gsrc(ord(i)).valid)>3
    istpos(:,n) = loc_osls([mic(:).pos],gsrc(ord(i)).tdoa*vs/Fs,gsrc(ord(i)).valid,gsrc(ord(i)).tdoaq);
    istqual(n) = gsrc(ord(i)).qual;
    n=n+1;
  end
end

[istpos;istqual]
