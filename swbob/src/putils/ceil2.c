int ceil2(int *n)
{
  int p=1;
  while(p<*n)
    p=p*2;
  return(p);
}

int ceil2_(int *n){
    int p;
    p=ceil2(n);
    return(p);
    }

