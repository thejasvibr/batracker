/*
  FUNCTION [idx] = ctrip(ref,fpair,fidx)
    She, 26.08.05

  Suche 2 Tripel aus f bei denen je ein anderes Paar gleiches Maximum liefert
  wie ref und die zueinander ein weiteres passendes Paar haben
  und fasse anschliessend Kombinationen zusammen

  Bsp: ref 1-2-3
       f1  1-2-4
       f2  1-3-4

  ref ist Vektor [1:MX] mit Indizes
  liefert Indizes des laengsten Satzes zurueck
  
  Achtung: Indizes beginnen mit 0 und sind vom Typ INT32
*/

#include <stdlib.h>
#include <mex.h>
#include <matrix.h>

unsigned int combtrip(const unsigned int *ref, unsigned int mx, const unsigned int *fpair, const unsigned int *fidx, unsigned int fanz, unsigned int *idx);
unsigned int merge(const unsigned int a[], const unsigned int b[], const unsigned int na, const unsigned int nb, unsigned int c[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* Schnittstelle zu MATLAB */
{
  unsigned int mx, fanz, i, n;
  unsigned int *ref, *fpair, *fidx;
  const double *mref, *mfpair, *mfidx;
  unsigned int *idx;

  /* Eingabewerte */
  mref = mxGetPr(prhs[0]);
  mx = mxGetNumberOfElements(prhs[0]);
  mfpair = mxGetPr(prhs[1]);
  mfidx = mxGetPr(prhs[2]);
  fanz = mxGetNumberOfElements(prhs[1])/3;
      
/* Array fuer Indizes */

  idx = mxCalloc(fanz,sizeof(int));

  ref = mxCalloc(mx,sizeof(int));
  fpair = mxCalloc(3*fanz,sizeof(int));
  fidx = mxCalloc(3*fanz,sizeof(int));

  for (i=0;i<mx;i++) 
    ref[i]=(int)mref[i]-1;
  for (i=0;i<3*fanz;i++) 
  {
    fpair[i]=(int)mfpair[i]-1;
    fidx[i]=(int)mfidx[i]-1;
  }

  n = combtrip(ref,mx,fpair,fidx,fanz,idx); 
  
  for (i=0;i<n;i++)
    idx[i]++;

  mxFree(ref);
  mxFree(fpair);
  mxFree(fidx);
    
/* Rueckgabewerte */

  plhs[0] = mxCreateNumericMatrix(1,n,mxINT32_CLASS,mxREAL);
  mxSetData(plhs[0],idx);
}

unsigned int combtrip(const unsigned int *ref, unsigned int mx, const unsigned int *fpair, const unsigned int *fidx, unsigned int fanz, unsigned int *idx)
{ 
/* 
  ref, idx: Vektoren der Laenge mx, 
            jede Stelle (0..mx-1) entspricht einem Mikrofonpaar,
            Inhalt ist jeweils die Nummer der untersuchten Laufzeitdifferenz (0,1,..)
            Wert -1 bedeutet, dass es (noch) keine passende Laufzeitdifferenz fuer dieses Paar gibt
            ref zu Beginn der Routine, idx am Ende
            
  fpair, fidx: Vektoren der Laenge 3*fanz
               fpair[3*n],fpair[3*n+1] und fpair[3*n+2] enthalten die Nummern von 3 Mikrofonpaaren (0..mx-1),
               die zusammen ein Tripel bilden. In fidx stehen die zugehoerigen Laufzeitdifferenznummern
            
  Variablennamen
  ct:    common tripel
  nct:   number of common tripels
  comb:  combinations
  ncomb: number of combinations
  fp:    free pairs
*/  
  
  unsigned int t1,t2,p1,p2,nct,ncomb,ct[fanz];
  unsigned int n,i,p;
  char act;
  struct fpt { unsigned int pair[3]; unsigned int idx[3]; unsigned int n;} fp[fanz];
  struct combt { unsigned int tripnr[fanz]; unsigned int n;} comb[fanz*(fanz-1)/2]; 
    
/* 
  Tripel aus f, die ein oder zwei gemeinsame Paare mit ref haben in fp speichern 
*/
  nct = 0;
  for (t1=0;t1<fanz;t1++)
  {
    if (ref[fpair[3*t1]]==fidx[3*t1] && ref[fpair[3*t1+1]]==fidx[3*t1+1] && ref[fpair[3*t1+2]]==-1)
    {
      ct[nct] = t1;
      fp[nct].pair[0] = fpair[3*t1+2];
      fp[nct].idx[0] = fidx[3*t1+2];
      fp[nct++].n = 1;
    }
    else if (ref[fpair[3*t1]]==fidx[3*t1] && ref[fpair[3*t1+1]]==-1 && ref[fpair[3*t1+2]]==fidx[3*t1+2])
    {
      ct[nct] = t1;
      fp[nct].pair[0] = fpair[3*t1+1];
      fp[nct].idx[0] = fidx[3*t1+1];
      fp[nct++].n = 1;
    }
    else if (ref[fpair[3*t1]]==-1 && ref[fpair[3*t1+1]]==fidx[3*t1+1] && ref[fpair[3*t1+2]]==fidx[3*t1+2])
    {
      ct[nct] = t1;
      fp[nct].pair[0] = fpair[3*t1];
      fp[nct].idx[0] = fidx[3*t1];
      fp[nct++].n = 1;
    }
    else if (ref[fpair[3*t1]]==fidx[3*t1] && ref[fpair[3*t1+1]]==-1 && ref[fpair[3*t1+2]]==-1)
    {
      ct[nct] = t1;
      fp[nct].pair[0] = fpair[3*t1+1];
      fp[nct].idx[0] = fidx[3*t1+1];
      fp[nct].pair[1] = fpair[3*t1+2];
      fp[nct].idx[1] = fidx[3*t1+2];
      fp[nct++].n = 2;
    }
    else if (ref[fpair[3*t1]]==-1 && ref[fpair[3*t1+1]]==fidx[3*t1+1] && ref[fpair[3*t1+2]]==-1)
    {
      ct[nct] = t1;
      fp[nct].pair[0] = fpair[3*t1];
      fp[nct].idx[0] = fidx[3*t1];
      fp[nct].pair[1] = fpair[3*t1+2];
      fp[nct].idx[1] = fidx[3*t1+2];
      fp[nct++].n = 2;
    }
    else if (ref[fpair[3*t1]]==-1 && ref[fpair[3*t1+1]]==-1 && ref[fpair[3*t1+2]]==fidx[3*t1+2])
    {
      ct[nct] = t1;
      fp[nct].pair[0] = fpair[3*t1];
      fp[nct].idx[0] = fidx[3*t1];
      fp[nct].pair[1] = fpair[3*t1+1];
      fp[nct].idx[1] = fidx[3*t1+1];
      fp[nct++].n = 2;
    }
  }
  
    
/*
  suche Kombinationen von 2 Tripeln aus fp, die zueinander ein weiteres gemeinsames Paar haben
  und speichere in comb die Tripelnummern
*/
  ncomb = 0;
  if (nct>1)
  {
    for (t1=0;t1<nct-1;t1++)
      for (t2=t1+1;t2<nct;t2++)
/*
  untersuche auf gemeinsames Paar von t1 und t2: 1-4 
*/
        for (p1=0;p1<fp[t1].n;p1++)
          for (p2=0;p2<fp[t2].n;p2++)
            if (fp[t1].pair[p1]==fp[t2].pair[p2])
              if (fp[t1].idx[p1]==fp[t2].idx[p2])
/*
  Paar liefert gleiche Laufzeitdifferenznummer 
*/
              {
                comb[ncomb].tripnr[0] = ct[t1];
                comb[ncomb].tripnr[1] = ct[t2];
                comb[ncomb++].n = 2;
                break;
                break;
              }
  }

  
/*
  fasse Kombinationen von Tripelnummern zusammen 
*/

  if (ncomb!=0)
  {
    act = 1;
    while (act)
    {
      act = 0;
      t1 = 0;
      while (t1<ncomb-1)
      {
        t2 = t1+1;
        while (t2<ncomb)
        {
          for (p1=0;p1<comb[t1].n;p1++)
            for (p2=0;p2<comb[t2].n;p2++)
              if (comb[t1].tripnr[p1]==comb[t2].tripnr[p2])
              {
/*
  Kombinationen passen zusammen, zusammensortieren 
*/
                comb[t1].n = merge(comb[t1].tripnr,comb[t2].tripnr,comb[t1].n,comb[t2].n,comb[t1].tripnr);
                comb[t2].n = 0;
                act = 1;
                break;
                break;              
              }
          t2++;
        }
        t1++;
      }
    }

/*
  groessten Kombinationserfolg suchen 
*/

    p2 = 0;
    for (t1=0;t1<ncomb;t1++)
      if (p2<comb[t1].n)
      {
        p2 = comb[t1].n;
        p1 = t1;
      }
  
    for (t1=0;t1<comb[p1].n;t1++)
      idx[t1] = comb[p1].tripnr[t1];
    
    return comb[p1].n;
  }
  else
    return 0;
}

unsigned int merge(const unsigned int a[], const unsigned int b[], const unsigned int na, const unsigned int nb, unsigned int c[])
/* 
  micht und sortiert die (sortierten) Teillisten a und b mit Laenge na bzw nb nach c
  c muss ausreichend gross dimensioniert sein 
  liefert Laenge von c zurueck 
*/
{
  unsigned int ai,bi,ci,nc;
  ai = 0;
  bi = 0;
  ci = 0;
  nc = na+nb;
  while (ci<nc)
  {
    if (ai<na && bi<nb)
    {
      if (a[ai]==b[bi]) /* a und b identisch */
      { 
        c[ci++] = a[ai++];
        bi++;
        nc--;
      }
      else if (a[ai]<b[bi])
        c[ci++] = a[ai++];
      else
        c[ci++] = b[bi++];
    }   
    else if (ai>=na) /* alle a einsortiert */
      c[ci++] = b[bi++];
    else /* alle b einsortiert */
      c[ci++] = a[ai++];
  }
  return nc;
  
}
