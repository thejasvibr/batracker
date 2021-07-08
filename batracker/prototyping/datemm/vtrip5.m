% FUNCTION list = vtrip(rcm12,rcm23,rcm13,params)
%   She, 27.07.05
%
% Suche gueltige Tripel aus den TDOA-Listen rcm12, rcm13 und rcm23
% liefert Tripel von Listenindizes und Treffermass zurueck
%
% params .PF

function list = vtrip(rcm12,rcm23,rcm13,params)

pfl = length(params.PF);
pfm = ceil(pfl/2);
th = min(params.PF);

% Summenmatrix aus rcm12 und rcm23
m13 = repmat(rcm12,size(rcm23'))+repmat(rcm23',size(rcm12));
lv = min([0;m13(:);rcm13(:)]);

% Musterfunktion fuer Summe aus rcm13 
%   Indizes [lv ... ], Maxima pfm rechts davon
p13m = zeros(max([0;m13(:);rcm13(:)])+pfl-lv,1);
p13i = zeros(max([0;m13(:);rcm13(:)])+pfl-lv,1);
for i=length(rcm13):-1:1
  p13m(rcm13(i)-lv+1:rcm13(i)-lv+pfl) = p13m(rcm13(i)-lv+1:rcm13(i)-lv+pfl)'+params.PF;
  p13i(rcm13(i)-lv+1:rcm13(i)-lv+pfl) = i;
end

% Suche Summen in Musterfunktion
testi = p13i(m13-lv+pfm);
testm = p13m(m13-lv+pfm);
match = find(testm>th);
[h1,h2] = ind2sub(size(m13),match);

list = [h1,h2,testi(match),testm(match)];

return;
