% FUNCTION micn = propag_trp(mic,src)
%  - akustisches Raummodell - She, 01.06.05
%
% Simulation der Schall-Ausbreitung von den Quellen src zu den Mikrofonen mic
% beruecksichtigt Laufzeit, Amplitudenabnahme mit 1/r und Pattern von
% Quellen und Mikrofonen
% ueberprueft Sichtbarkeit
%
% in:  mic, src, vs, Fs, wavlen
% out: micn.wav

function micn = propag_tr(mic,src,wall)

% Modell: micsig(t)=sum_(src)(srcsig(t-delay)/dist*a_mic*a_src)

global vs Fs wavlen dmsg;
micn = mic;

for m = 1:length(micn)
  if dmsg disp(['mic ',num2str(m)]); end;
% Initialisiere Mikrofonsignal
  micn(m).wav = zeros(wavlen,1);

  for s = 1:length(src)
% Sichtbarkeitspruefung:
%   Suche Schnittpunkt der Geraden mic-src mit den Waenden, der am nahesten 
%   bei mic liegt, die zugehoerige Wand muss die letzte Wand in der Quellbezeichnung
%   sein. Verlege dann mic an den Reflexionspunkt und src an zu dieser Spiegelung
%   gehoerende Urbildquelle und pruefe erneut, solange bis Urbildquelle echte Quelle
%   ist.
    visible = true;
    wi = strfind(src(s).name,'_');
    moi = mic(m).pos;
    soi = s;
    while ~isempty(wi)
% Schnittpunkt der Geraden moi-src mit den Waenden, der am nahesten bei moi liegt
      lmin = 1;
      for w = 1:length(wall)
        ne = dot(src(soi).pos-moi,wall(w).nv);
        if abs(ne)>3*eps
          l = dot(wall(w).pos1-moi,wall(w).nv)/ne;
          if (l>0 & l<lmin) lmin = l; wmin = w; end
        end
      end
% Vergleich der Wandnamen
      if strcmp(src(soi).name(wi(end)+1:end),wall(wmin).name)
        moi = moi+lmin*(src(soi).pos-moi);
        for s2 = 1:length(src)
          if strcmp(src(s2).name,src(soi).name(1:wi(end)-1)) soi = s2; end
        end
        wi = wi(1:end-1);
      else
        visible = false;
        wi=[];
      end
    end
    
    if visible
      if dmsg>=2 disp(['Quelle ',src(s).name,' von ',mic(m).name,' sichtbar']); end;
% Abstand des Quell-Mikrofon-Paares
      dist = norm(micn(m).pos-src(s).pos);
      delay = round(dist/vs*Fs);
% Verzoegerung
      w = [zeros(delay,1);src(s).wav;zeros(wavlen,1)];
% Amplitudenabnahme aufgrund Ausbreitung und Verkuerzung
      w = 1/dist*w([1:wavlen],1);
% Amplitudenabnahme durch Richtcharakteristik der Quelle
      if isfield(src(s),'dir')
        d = mic(m).pos-src(s).pos;
        xi = acos(dot(src(s).dir,d)/sqrt(src(s).dir'*src(s).dir*d'*d));
        w = w*eval(src(s).eff);
      end
% Amplitudenabnahme durch Richtcharakteristik des Mikrofons
      if isfield(mic(m),'dir')
        d = src(s).pos-mic(m).pos;
        xi = acos(dot(mic(m).dir,d)/sqrt(mic(m).dir'*mic(m).dir*d'*d));
        w = w*eval(mic(m).eff);
      end
% Ueberlagerung
      micn(m).wav = micn(m).wav+w;
    end
  end
  maxamp(m) = max(micn(m).wav);
end
maxampt = max(maxamp);

% Normierung
if maxampt > 1 
  for m = 1:length(micn)
    micn(m).wav = micn(m).wav/maxampt;
  end
end

return;
