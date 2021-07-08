% FUNCTION est=osls(mic,delta,valid,qual)
%  05.05.06 She
%
% 3D-paarweises Lokalisieren einer Quelle mit n Sensoren
% one-step-least-squares nach Huang 
% aequivalent zur Kreisinterpolation
%
% mic:   (3xn)-Matrix mit Koordinaten der Sensoren
% delta: (nxn)-Matrix mit Abstandsdifferenzen
%
% angepasst auf Ausgabedaten von DATEMM

function est = osls(mic,delta,valid,qual)

i=1;
for m1=1:length(mic)-1
  for m2=m1+1:length(mic)
    pdelta(m1,m2)=delta(i);
    pdelta(m2,m1)=-delta(i);
    pvalid(m1,m2)=valid(i);
    pvalid(m2,m1)=valid(i);
    pqual(m1,m2)=qual(i);
    pqual(m2,m1)=qual(i);
    i=i+1;
  end
end

% Referenzsensor suchen und in Koordinatenursprung legen
dummy=find(sum(pvalid)==max(sum(pvalid)));
dummy2=find(sum(pqual(:,dummy))==max(sum(pqual(:,dummy))));
ref = dummy(dummy2(1));
if (sum(pvalid(:,ref))<3)
  warning('zu wenig TDOAs');
  est = [-Inf,-Inf,-Inf]';
  return;
end

nref = find(pvalid(:,ref)==1);
      
% Koordinaten der anderen Sensoren und Abstaende
S = (mic(:,nref)-mic(:,ref)*ones(1,length(nref)))';
d = pdelta(ref,nref)';

A = [S,d];

% Kleinste Quadrate-Schaetzung
t = .5*inv(A'*A)*A'*(S(:,1).^2+S(:,2).^2+S(:,3).^2-d.^2);
est = mic(:,ref)+t(1:end-1);
        
return
