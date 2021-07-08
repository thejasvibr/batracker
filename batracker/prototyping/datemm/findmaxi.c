/*
  FUNCTION idx = findmaxi(val,n,thu,thl)
    She, 29.07.05

  Suche Position der ersten Maxima in val,
  zunaechst alle, die groesser als thu*globalmax sind
  falls das weniger als n sind, dann weitere 
    die groeser als thl*globalmax sind, in Summe maximal n
    
  Achtung: Indizes beginnen mit 0 und sind vom Typ INT32
  Aufruf der Werte in MATLAB: val(double(findmaxi(val,n,thu,thl))+1)
*/

#include <stdlib.h>
#include <mex.h>
#include <matrix.h>

int findrelmax(const double *val, int *idx, int len, int n, double thu, double thl);
double globmax(const double *val, int len);
void quicksort(const double *val, int *idx, int l, int r);
int partition(const double *val, int *idx, int left, int right, double thres);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* Schnittstelle zu MATLAB */
{
  int n, len, i;
  const double *val;
  double thu, thl;
  int *idx;
  
/* Eingabewerte */
  val = mxGetPr(prhs[0]);
  len = mxGetNumberOfElements(prhs[0]);
  n = mxGetScalar(prhs[1]);
  thu = mxGetScalar(prhs[2]);
  thl = mxGetScalar(prhs[3]);
  
/* Array fuer Indizes */
  idx = mxCalloc(len,sizeof(int));
  
  i = findrelmax(val,idx,len,n,thu,thl);
  
/* Rueckgabewerte */
  plhs[0] = mxCreateNumericMatrix(i,1,mxINT32_CLASS,mxREAL);
  mxSetData(plhs[0],idx);
}


int findrelmax(const double *val, int *idx, int len, int n, double thu, double thl)
/* Hauptroutine */
/* Verfahren siehe oben */
{
  int i, nthu, nthl;
  double max, dummy;
  
/* initialisiere Indizes-Array */
  for (i=0; i<len; i++) idx[i] = i;
  
  max = globmax(val,len);
  
/* spalte Werte > thu*globalmax ab und sortiere diese */
  nthu = partition(val,idx,0,len-1,thu*max);   
  
  dummy = idx[nthu]; idx[nthu] = idx[len-1]; idx[len-1] = dummy;
  quicksort(val,idx,0,nthu);
  
  if (nthu < n)
  {
    nthl = partition(val,idx,0,len-1,thl*max);   
    dummy = idx[nthl]; idx[nthl] = idx[len-1]; idx[len-1] = dummy;

    quicksort(val,idx,nthu,nthl);
    if (nthl < n)
      return nthl;
    else
      return n;
  }
  else 
    return n;
}

int partition(const double *val, int *idx, int left, int right, double thres)
/* Quicksort-Partitionierung 
     sortiere alle Werte von val, die groesser als thres sind, auf die linke
     Seite von thres und alle kleineren auf die rechte Seite. Liefere die
     Position des Werts thres zurueck. 
     
   insofern abgewandelt, als dass val unveraendert bleibt und nur die Indizierung
   idx veraendert wird. idx wird nur im Bereich left..right betrachtet. */
{
  int dummy, i;
  left--;
  for (;;)
  {
    while ((val[idx[++left]] > thres)&&(left<right));
    while ((val[idx[--right]] < thres)&&(left<right));
    if (left >= right) break;
    dummy = idx[left]; idx[left] = idx[right]; idx[right] = dummy;
  }
  return left;
}


void quicksort(const double *val, int *idx, int left, int right)
/* Quicksort: rekursive Sortierung von val im Bereich idx=left..right */
{
  int pivot, dummy;
  if (left < right)
  {
    pivot = partition(val,idx,left,right,val[idx[right]]);   
    dummy = idx[pivot]; idx[pivot] = idx[right]; idx[right] = dummy;
    quicksort(val,idx,left,pivot-1);
    quicksort(val,idx,pivot+1,right);
  }
}

double globmax(const double *val, int len)
/* liefert maximalen Wert des Arrays val mit Laenge len */  
{
  int i;
  double m;
  m = val[0];
  for (i=1; i<len; i++)
    if (val[i]>m) m=val[i];
  return m;  
}

