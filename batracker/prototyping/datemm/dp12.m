% FUNCTION [idx,val] = dp(xc21,acm1,acmv1,acm2,acmv2,mdist,params)
%   She, 31.08.05
%
% Suche Direktpfadlaufzeiten in Kreuzkorrelierter xc21
% mit Hilfe der Autokorrelierten ac1 und ac2
% durch Abstandssuche
% mit Aufwertung der direkten Pfade und Abwertung der Echopfade nach Autokorrelationsbetrag
% beruecksichtigt physikalisch maximal moeglichen Abstandsbetrag mdist 
% naheliegende Maxima (Abstand <= .DIST) werden zusammengefasst.
%
% Indizes Autokorrelierte [0 1 2 ... N-1]
% Indizes Kreuzkorrelierte [-N+1 ... -1 0 1 2 ... N-1]
% params .PF, .AUF, .AB, .FIND, .NX, .TUX, .TLA, .DIST
%

function [idx,val] = dp(xc21,acm1,acmv1,acm2,acmv2,mdist,params)
  
global DMSG;
    
N = .5*(length(xc21)+1);
pfm = length(params.PF)-1;
  
% relevante Betragsmaxima der Kreuzkorrelierten
xcm = double(feval(params.FIND,abs(xc21),params.NX,params.TUX,params.TLX))+1;
xcmv = abs(xc21(xcm));
dummy = sortrows([xcm,xcmv]);

% benachbarte Maxima zusammenfassen
xcm = [];
xcmv = [];
i = 1;
n1 = 1;
while n1<=size(dummy,1)
  n2 = 1;
  while (n1+n2<=size(dummy,1)) & (dummy(n1+n2,1)-dummy(n1+n2-1,1)<=params.DIST)
    n2 = n2+1;
  end
  xcm(i,1) = ceil(mean(dummy(n1:n1+n2-1,1)));
  xcmv(i,1) = sum(dummy(n1:n1+n2-1,2));
  i = i+1;
  n1 = n1+n2;
end

% passende Direkt-/Echopfade
mdp = zeros(size(xcm));
mep = zeros(size(xcm));

% Korrelationsmaxima-Differenzen
xcfp = 0;
for n1=1:size(xcm)-1
  for n2=n1+1:size(xcm)
    xcfp = xcfp+1;
    xcf(xcfp) = xcm(n2)-xcm(n1);
    xcfn1(xcfp) = n1;
    xcfn2(xcfp) = n2;
  end
end

% Rastersuche
for i=1:length(acm1)
  for n=1:xcfp
    if abs(xcf(n)-acm1(i))<=pfm
      mdp(xcfn1(n)) = mdp(xcfn1(n))+acmv1(i);
      mep(xcfn2(n)) = mep(xcfn2(n))+acmv1(i).*params.PF(abs(xcf(n)-acm1(i))+1)';
    end
  end
end
for i=1:length(acm2)
  for n=1:xcfp
    if abs(xcf(n)-acm2(i))<=pfm
      mdp(xcfn2(n)) = mdp(xcfn2(n))+acmv2(i);
      mep(xcfn1(n)) = mep(xcfn1(n))+acmv2(i).*params.PF(abs(xcf(n)-acm2(i))+1)';
    end
  end
end

% Auswirkung auf Korrelationsmaxima
xcmi = find(xcmv+params.AUF*mdp+params.AB*mep>min(xcmv));
dummy = xcm(xcmi(find(abs(xcm(xcmi)-N)<mdist)));


idx = dummy-N;
val = xc21(dummy);

return;
