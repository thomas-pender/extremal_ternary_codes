# include "./search8.h"

typedef uint_fast32_t UINT;

// misc //////////////////////////////////////////////////////////////////////////

void err_display () __attribute__ ((noreturn));
void err_display(int throw,
                 char const name[static 1],
                 char const file[static 1],
                 char const func[static 1],
                 unsigned long line)
{
  fprintf(stderr,"ERROR -- %d: %s by %s\n\t%s:%s:%lu\n",
          throw,strerror(throw),name,file,func,line);
  exit(EXIT_FAILURE);
}

void err_short () __attribute__ ((noreturn));
void err_short(char const name[static 1],
               char const file[static 1],
               char const func[static 1],
               unsigned long line)
{
  fprintf(stderr,"ERROR -- by %s\n\t%s:%s:%lu\n",
          name,file,func,line);
  exit(EXIT_FAILURE);
}

static inline
void swap_num(unsigned a[restrict static 1], unsigned b[restrict static 1])
{
  unsigned temp=*a;
  *a=*b;
  *b=temp;
}

static inline
void rev_array(unsigned a[static 1], size_t left, size_t right)
{
  for(; left<right; ++left, --right)
    swap_num(&a[left],&a[right]);
}

static inline
unsigned min(unsigned a, unsigned b)
{
  if(a <= b) return a;
  else return b;
}

static inline
void copy_mat(UINT *restrict A[restrict static 1],
              UINT *restrict B[restrict static 1])
{
  size_t i,j,len=8*LEN,Len=16*LEN;
  for(i=0; i<len; ++i)
    for(j=0; j<Len; B[i][j]=A[i][j],++j);
}

void print_mat(UINT **A, size_t nrows, size_t ncols)
{
  size_t i,j;
  for(i=0; i<nrows; ++i){
    for(j=0; j<ncols; printf("%2lu",A[i][j++]));
    printf("\n");
  }
}

static inline
void transpose(UINT *restrict A[restrict static 1],
               UINT *restrict B[restrict static 1])
{
  size_t i,j,len=8*LEN,Len=16*LEN;
  for(i=0; i<Len; ++i)
    for(j=0; j<len; B[i][j]=A[j][i],++j);
}

// sequence arithmetic ///////////////////////////////////////////////////////////

void corrs(UINT a[restrict static 1],
           UINT corr[restrict static 1],
           unsigned clen)
{
  size_t i,j;
  UINT sum;
  for(j=1; j<=clen; ++j){
    sum=0;
# if NEGA
    for(i=0; i<LEN-j; sum+=a[i]*a[i+j],++i);
    for(i=0; i<j; sum+=(2*a[i]*a[LEN+i-j])%3,++i);
# else
    for(i=0; i<LEN; sum+=a[i]*a[(i+j)%LEN],++i);
# endif
    corr[j-1]=sum%3;
  }
}

UINT *circ_mult(UINT a[restrict static 1], UINT b[restrict static 1])
{
  size_t i,j;
  UINT sum,*mult=(UINT*)malloc(LEN*sizeof(UINT));
  for(i=0,sum=0; i<LEN; sum+=a[i]*b[i],++i);
  mult[0]=sum%3;
# if NEGA
  for(j=LEN-1; j; --j){
    sum=0;
    for(i=0; i<LEN-j; sum+=2*a[i]*b[i+j],++i);
    for(i=LEN-j; i<LEN; sum+=a[i]*b[(i+j)%LEN],++i);
    mult[LEN-j]=sum%3;
  }
# else
  for(j=LEN-1; j; --j){
    sum=0;
    for(i=0; i<LEN; sum+=a[i]*b[(i+j)%LEN],++i);
    mult[LEN-j]=sum%3;
  }
# endif
  return mult;
}

UINT *circ_add(UINT a[restrict static 1],
               UINT b[restrict static 1],
               UINT c[restrict static 1],
               UINT d[restrict static 1])
{
  size_t i;
  UINT *mult=(UINT*)malloc(LEN*sizeof(UINT));
  for(i=0; i<LEN; ++i)
    mult[i]=(a[i]+b[i]+c[i]+d[i])%3;
  return mult;
}

UINT *circ_minus(UINT a[restrict static 1],
                 UINT b[restrict static 1])
{
  size_t i;
  UINT *sub=(UINT*)malloc(LEN*sizeof(UINT));
  for(i=0; i<LEN; ++i)
    sub[i]=(a[i] + 2*b[i])%3;
  return sub;
}

// generator matrix //////////////////////////////////////////////////////////////

UINT **circ(UINT a[restrict static 1], unsigned transpose, unsigned minus)
{
  size_t i,j;
  UINT **A=(UINT**)malloc(LEN*sizeof(UINT*));
  for(i=0; i<LEN; A[i++]=(UINT*)malloc(LEN*sizeof(UINT)));
  if(transpose && minus){
    A[0][0]=(2*a[0])%3;
# if NEGA
    for(i=1; i<LEN; A[0][i]=a[LEN-i],++i);
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j) A[i][j]=A[i-1][j-1];
        else A[i][j]=(2*A[i-1][LEN-1])%3;
      }
# else
    for(i=1; i<LEN; A[0][i]=(2*a[LEN-i])%3,++i);
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j) A[i][j]=A[i-1][j-1];
        else A[i][j]=A[i-1][LEN-1];
      }
# endif
  }
  else if(transpose && !minus){
    A[0][0]=a[0];
# if NEGA
    for(i=1; i<LEN; A[0][i]=(2*a[LEN-i])%3,++i);
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j) A[i][j]=A[i-1][j-1];
        else A[i][j]=(2*A[i-1][LEN-1])%3;
      }
# else
    for(i=1; i<LEN; A[0][i]=a[LEN-i],++i);
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j) A[i][j]=A[i-1][j-1];
        else A[i][j]=A[i-1][LEN-1];
      }
# endif
  }
  else if(!transpose && minus){
    for(i=0; i<LEN; A[0][i]=(2*a[i])%3,++i);
# if NEGA
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j) A[i][j]=A[i-1][j-1];
        else A[i][j]=(2*A[i-1][LEN-1])%3;
      }
# else
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j) A[i][j]=A[i-1][j-1];
        else A[i][j]=A[i-1][LEN-1];
      }
# endif
  }
  else{
    for(i=0; i<LEN; A[0][i]=a[i],++i);
# if NEGA
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j) A[i][j]=A[i-1][j-1];
        else A[i][j]=(2*A[i-1][LEN-1])%3;
      }
# else
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j) A[i][j]=A[i-1][j-1];
        else A[i][j]=A[i-1][LEN-1];
      }
# endif
  }
  return A;
}

UINT **bcirc(UINT a[restrict static 1], unsigned transpose, unsigned minus)
{
  size_t i,j;
  UINT **A=(UINT**)malloc(LEN*sizeof(UINT*));
  for(i=0; i<LEN; A[i++]=(UINT*)malloc(LEN*sizeof(UINT)));
  if(transpose && minus){
    A[0][LEN-1]=(2*a[0])%3;
# if NEGA
    for(i=0; i<LEN-1; A[0][i]=a[i+1],++i);
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j == LEN-1) A[i][j]=(2*A[i-1][0])%3;
        else A[i][j]=A[i-1][j+1];
      }
# else
    for(i=0; i<LEN-1; A[0][i]=(2*a[i+1])%3,++i);
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j == LEN-1) A[i][j]=A[i-1][0];
        else A[i][j]=A[i-1][j+1];
      }
# endif
  }
  else if(transpose && !minus){
    A[0][LEN-1]=a[0];
# if NEGA
    for(i=0; i<LEN-1; A[0][i]=(2*a[i+1])%3,++i);
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j == LEN-1) A[i][j]=(2*A[i-1][0])%3;
        else A[i][j]=A[i-1][j+1];
      }
# else
    for(i=0; i<LEN-1; A[0][i]=a[i+1],++i);
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j == LEN-1) A[i][j]=A[i-1][0];
        else A[i][j]=A[i-1][j+1];
      }
# endif
  }
  else if(!transpose && minus){
    for(i=0; i<LEN; A[0][i]=(2*a[LEN-i-1])%3,++i);
# if NEGA
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j == LEN-1) A[i][j]=(2*A[i-1][0])%3;
        else A[i][j]=A[i-1][j+1];
      }
# else
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j == LEN-1) A[i][j]=A[i-1][0];
        else A[i][j]=A[i-1][j+1];
      }
# endif
  }
  else{
    for(i=0; i<LEN; A[0][i]=a[LEN-i-1],++i);
# if NEGA
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j == LEN-1) A[i][j]=(2*A[i-1][0])%3;
        else A[i][j]=A[i-1][j+1];
      }
# else
    for(i=1; i<LEN; ++i)
      for(j=0; j<LEN; ++j){
        if(j == LEN-1) A[i][j]=A[i-1][0];
        else A[i][j]=A[i-1][j+1];
      }
# endif
  }
  return A;
}

UINT **gen_mat(UINT a1[restrict static 1],
               UINT a2[restrict static 1],
               UINT a3[restrict static 1],
               UINT a4[restrict static 1],
               UINT a5[restrict static 1],
               UINT a6[restrict static 1],
               UINT a7[restrict static 1],
               UINT a8[restrict static 1])
{
  size_t i,j,len=8*LEN,Len=16*LEN;
  UINT **A1=circ(a1,0,0),
    **A2=circ(a2,0,0),**MA2=circ(a2,0,1),
    **A3=bcirc(a3,0,0),**MA3=bcirc(a3,0,1),**A3T=bcirc(a3,1,0),**MA3T=bcirc(a3,1,1),
    **A4=bcirc(a4,0,0),**MA4=bcirc(a4,0,1),**A4T=bcirc(a4,1,0),**MA4T=bcirc(a4,1,1),
    **A5=bcirc(a5,0,0),**MA5=bcirc(a5,0,1),**A5T=bcirc(a5,1,0),**MA5T=bcirc(a5,1,1),
    **A6=bcirc(a6,0,0),**MA6=bcirc(a6,0,1),**A6T=bcirc(a6,1,0),**MA6T=bcirc(a6,1,1),
    **A7=bcirc(a7,0,0),**MA7=bcirc(a7,0,1),**A7T=bcirc(a7,1,0),**MA7T=bcirc(a7,1,1),
    **A8=bcirc(a8,0,0),**MA8=bcirc(a8,0,1),**A8T=bcirc(a8,1,0),**MA8T=bcirc(a8,1,1),
    **P[DIM][DIM]={{A1,A2,A4,A3,A6,A5,A8,A7,},
                   {MA2,A1,A3,MA4,A5,MA6,A7,MA8,},
                   {MA4,MA3,A1,A2,MA8T,A7T,A6T,MA5T,},
                   {MA3,A4,MA2,A1,A7T,A8T,MA5T,MA6T,},
                   {MA6,MA5,A8T,MA7T,A1,A2,MA4T,A3T,},
                   {MA5,A6,MA7T,MA8T,MA2,A1,A3T,A4T,},
                   {MA8,MA7,MA6T,A5T,A4T,MA3T,A1,A2,},
                   {MA7,A8,A5T,A6T,MA3T,MA4T,MA2,A1,},},
    **G=(UINT**)malloc(len*sizeof(UINT*));
  for(i=0; i<len; G[i++]=(UINT*)malloc(Len*sizeof(UINT)));

  for(i=0; i<len; ++i)
    for(j=0; j<len; ++j){
      if(i == j) G[i][j]=1;
      else G[i][j]=0;
    }
  for(i=0; i<len; ++i){
    size_t q=(i-(i%LEN))/LEN,r=i%LEN;
    for(j=0; j<len; ++j)
      G[i][len+j]=P[q][(j-(j%LEN))/LEN][r][j%LEN];
  }

  for(i=0; i<LEN; ++i){
    free(A1[i]); free(A2[i]); free(MA2[i]);
    free(A3[i]); free(MA3[i]); free(A3T[i]); free(MA3T[i]);
    free(A4[i]); free(MA4[i]); free(A4T[i]); free(MA4T[i]);
    free(A5[i]); free(MA5[i]); free(A5T[i]); free(MA5T[i]);
    free(A6[i]); free(MA6[i]); free(A6T[i]); free(MA6T[i]);
    free(A7[i]); free(MA7[i]); free(A7T[i]); free(MA7T[i]);
    free(A8[i]); free(MA8[i]); free(A8T[i]); free(MA8T[i]);
  }
  free(A1); free(A2); free(MA2);
  free(A3); free(MA3); free(A3T); free(MA3T);
  free(A4); free(MA4); free(A4T); free(MA4T);
  free(A5); free(MA5); free(A5T); free(MA5T);
  free(A6); free(MA6); free(A6T); free(MA6T);
  free(A7); free(MA7); free(A7T); free(MA7T);
  free(A8); free(MA8); free(A8T); free(MA8T);

  return G;
}

// rref //////////////////////////////////////////////////////////////////////////

static inline
void swap_rows(UINT **A, size_t r1, size_t r2)
{
  UINT *temp=A[r1];
  A[r1]=A[r2];
  A[r2]=temp;
}

static inline
void swap_cols(UINT **A, size_t nrows, size_t c1, size_t c2)
{
  UINT temp;
  for(size_t i=0; i<nrows; ++i){
    temp=A[i][c1];
    A[i][c1]=A[i][c2];
    A[i][c2]=A[i][c1];
  }
}

static inline
void row_mult(UINT **A, size_t ncols, size_t r, UINT mult)
{
  for(size_t i=0; i<ncols; A[r][i]=(mult*A[r][i])%3,++i);
}

static inline
unsigned row_mult_add(UINT **A, size_t ncols, size_t r1, size_t r2, UINT mult)
{
  unsigned check=1,len=8*LEN;
  for(size_t i=0; i<ncols; ++i){
    A[r2][i]=(A[r2][i] + (mult*A[r1][i]))%3;
    if(i>=len && check && A[r2][i]) --check;
  }
  return check;
}

UINT **rref(UINT *restrict G[restrict static 1])
{
  size_t i,j,k,nrows=8*LEN,ncols=16*LEN,row,col;
  unsigned throw;
  UINT **A=(UINT**)malloc(nrows*sizeof(UINT*)),mult;
  for(i=0; i<nrows; A[i++]=(UINT*)malloc(ncols*sizeof(UINT)));
  copy_mat(G,A);

  for(i=0,col=nrows,row=nrows; i<row; ++i,++col){
    if(col == ncols){
      row=i;
      break;
    }

    if(!A[i][col]){
      for(j=i+1; j<row; ++j)
        if(A[j][col]) break;

      if(j == row){
        for(k=col+1; k<ncols; ++k)
          if(A[i][k]) break;

        if(k == ncols){
          swap_rows(A,i,row-1);
          --row;
          --i;
          --col;
        }
        else
          swap_cols(A,row,col,k);
      }
      else
        swap_rows(A,i,j);
    }

    if(A[i][col] == 2) row_mult(A,ncols,i,2);

    for(j=i+1; j<row; ++j)
      if(A[j][col]){
        mult=3-A[j][col];
        if(row_mult_add(A,ncols,i,j,mult)){
          swap_rows(A,j,row-1);
          --row;
          --j;
        }
      }
  }

  for(i=1,col=nrows+1; i<row && col<ncols; ++i,++col)
    for(j=0; j<i; ++j)
      if(A[j][col]){
        mult=3-A[j][col];
        throw=row_mult_add(A,ncols,i,j,mult);
      }

  return A;
}

// find min weight ///////////////////////////////////////////////////////////////

unsigned binsearch(UINT a[static 1], size_t left, size_t right, UINT num)
{
  size_t mid;
  while(left<=right){
    mid=left+(right-left)/2;
    if(a[mid] == num) return 1U;
    if(a[mid] < num) left=mid+1;
    else right=mid-1;
  }
  return 0U;
}

unsigned weight(UINT a[restrict static 1], UINT *restrict A[restrict static 1])
{
  size_t i,j,len=8*LEN,Len=16*LEN;
  UINT sum;
  unsigned wt=0;
  for(i=0; i<Len; ++i){
    sum=0;
    for(j=0; j<len; sum+=a[j]*A[i][j],++j);
    if(sum%3) ++wt;
  }
  return wt;
}

void vectors(unsigned k,
             UINT *restrict A[restrict static 1],
             UINT *restrict B[restrict static 1],
             unsigned wt[static 1])
{
  size_t i,j,len=8*LEN;
  unsigned w1=len+1,w2=len+1;
  UINT counter,lim=(((UINT)1)<<k),
    dat[k],v[len],num;
  for(i=0; i<len; v[i++]=0);

  for(i=0; i<k; dat[i]=(UINT)i,++i);
  for(counter=0; counter<lim; ++counter){
    for(i=0; i<k; ++i){
      if(counter & (((UINT)1)<<i)) v[dat[i]]=2;
      else v[dat[i]]=1;
    }
    w1=min(w1,weight(v,A));
    w2=min(w2,weight(v,B));
  }
  for(i=0; i<len; v[i++]=0);

  for(unsigned check=1; check;){
    for(i=0; i<k; ++i){
      num=dat[i]+1;
      if(!binsearch(dat,i+1,k-1,num)){
        if(i==k-1 && num==len){
          --check;
          break;
        }
        ++dat[i];
        for(j=0; j<i; dat[j]=(UINT)j,++j);
        for(counter=0; counter<lim; ++counter){
          for(j=0; j<k; ++j){
            if(counter & (((UINT)1)<<j)) v[dat[j]]=2;
            else v[dat[j]]=1;
          }
          w1=min(w1,weight(v,A));
          w2=min(w2,weight(v,B));
        }
        for(j=0; j<len; v[j++]=0);
        break;
      }
    }
  }

  *wt=min(*wt,w1);
  *wt=min(*wt,w2);
}

unsigned min_wt(UINT *restrict A[restrict static 1])
{
  size_t i,len=8*LEN,Len=16*LEN;
  UINT **B=rref(A),
    **AT=(UINT**)malloc(Len*sizeof(UINT*)),
    **BT=(UINT**)malloc(Len*sizeof(UINT*));
  for(i=0; i<Len; ++i){
    AT[i]=(UINT*)malloc(len*sizeof(UINT));
    BT[i]=(UINT*)malloc(len*sizeof(UINT));
  }
  transpose(A,AT);
  transpose(B,BT);
  for(i=0; i<len; free(B[i]),++i);
  free(B);

  unsigned lb=1,ub=len+1;
  for(unsigned w=1; lb<ub; ++w){
    vectors(w,AT,BT,&ub);
    lb=2*(w+1);
  }

  for(i=0; i<Len; ++i){
    free(AT[i]);
    free(BT[i]);
  }
  free(AT);
  free(BT);

  return ub;
}

// amicability check /////////////////////////////////////////////////////////////

void print_array(UINT a[static 1])
{
  printf("(");
  for(size_t i=0; i<LEN; ++i)
    printf("%-2lu",a[i]);
  printf("\b)\n");
}

unsigned check_amicable(UINT a1[restrict static 1],
                        UINT a2[restrict static 1],
                        UINT a3[restrict static 1],
                        UINT a4[restrict static 1],
                        UINT a5[restrict static 1],
                        UINT a6[restrict static 1],
                        UINT a7[restrict static 1],
                        UINT a8[restrict static 1])
{
  size_t i;
  UINT *A2A1T=circ_mult(a2,a1),
    *A1A2T=circ_mult(a1,a2),
    *A4A3T=circ_mult(a4,a3),
    *A3A4T=circ_mult(a3,a4),
    *A6A5T=circ_mult(a6,a5),
    *A5A6T=circ_mult(a5,a6),
    *A8A7T=circ_mult(a8,a7),
    *A7A8T=circ_mult(a7,a8),
    *T1=circ_minus(A2A1T,A1A2T),
    *T2=circ_minus(A4A3T,A3A4T),
    *T3=circ_minus(A6A5T,A5A6T),
    *T4=circ_minus(A8A7T,A7A8T),
    *T=circ_add(T1,T2,T3,T4);
  unsigned check;

  for(i=0,check=1; i<LEN; ++i)
    if(T[i]){
      --check;
      break;
    }

  free(A2A1T); free(A1A2T);
  free(A4A3T); free(A3A4T);
  free(A6A5T); free(A5A6T);
  free(A8A7T); free(A7A8T);
  free(T1); free(T2);
  free(T3); free(T4);
  free(T);

  return check;
}


UINT** permutations(UINT a1[restrict static 1],
                    UINT a2[restrict static 1],
                    UINT a3[restrict static 1],
                    UINT a4[restrict static 1],
                    UINT a5[restrict static 1],
                    UINT a6[restrict static 1],
                    UINT a7[restrict static 1],
                    UINT a8[restrict static 1])
{
  size_t i,j;
  unsigned dat[DIM];
  UINT *P[DIM]={a1,a2,a3,a4,a5,a6,a7,a8,};

  for(i=0; i<DIM; ++i)
    dat[i]=(unsigned)i;
  if(check_amicable(P[0],P[1],P[2],P[3],P[4],P[5],P[6],P[7]))
    return gen_mat(P[0],P[1],P[2],P[3],P[4],P[5],P[6],P[7]);

  for(;;){
    for(i=DIM-1; i; --i)
      if(dat[i-1] < dat[i]) break;
    if(!i) break;
    else{
      for(j=DIM-1; j>i; --j)
        if(dat[j] > dat[i-1]) break;
      swap_num(&dat[j],&dat[i-1]);
      rev_array(dat,i,DIM-1);
      if(check_amicable(P[dat[0]],P[dat[1]],P[dat[2]],P[dat[3]],P[dat[4]],P[dat[5]],P[dat[6]],P[dat[7]]))
        return gen_mat(P[dat[0]],P[dat[1]],P[dat[2]],P[dat[3]],P[dat[4]],P[dat[5]],P[dat[6]],P[dat[7]]);
    }
  }

  return NULL;
}

// partitions ////////////////////////////////////////////////////////////////////

void partitions(int n,
                UINT *restrict*restrict sub[restrict static 1],
                UINT *restrict*restrict corr[restrict static 1],
                UINT lengths[restrict static 1],
                unsigned clen)
{
  register size_t i,a,b,c,d,e,f,g,h,len=8*LEN,Len=16*LEN;
  unsigned check,D=3*floor(4*LEN/3.0)+3;
  UINT sum,**G,
    *z=(UINT*)calloc(LEN,sizeof(UINT)),
    *zz=(UINT*)calloc(LEN,sizeof(UINT)),
    *zzz=(UINT*)calloc(LEN,sizeof(UINT));

  int p[n];
  int k=0;
  p[0]=n;

  while(1){
    if(k==4 && p[0]<=LEN){
      for(a=0; a<lengths[p[0]-1]; ++a)
        for(b=0; b<lengths[p[1]-1]; ++b)
          for(c=0; c<lengths[p[2]-1]; ++c)
            for(d=0; d<lengths[p[3]-1]; ++d)
              for(e=0; e<lengths[p[4]-1]; ++e){
                check=1;
                for(i=0; i<clen; ++i){
                  sum=corr[p[0]-1][a][i]+
                    corr[p[1]-1][b][i]+
                    corr[p[2]-1][c][i]+
                    corr[p[3]-1][d][i]+
                    corr[p[4]-1][e][i];
                  if(sum%3){
                    --check;
                    break;
                  }
                }
                if(check){
                  G=permutations(z,zz,zzz,sub[p[0]-1][a],sub[p[1]-1][b],sub[p[2]-1][c],sub[p[3]-1][d],sub[p[4]-1][e]);
                  if(G!=NULL && D==min_wt(G)){
                    for(i=len; i<Len; ++i)
                      printf("%2lu",G[0][i]);
                    printf("\n");
                  }
                  if(G!=NULL){
                    for(i=0; i<len; ++i)
                      free(G[i]);
                    free(G);
                  }
                }
              }
    }
    if(k==5 && p[0]<=LEN){
      for(a=0; a<lengths[p[0]-1]; ++a)
        for(b=0; b<lengths[p[1]-1]; ++b)
          for(c=0; c<lengths[p[2]-1]; ++c)
            for(d=0; d<lengths[p[3]-1]; ++d)
              for(e=0; e<lengths[p[4]-1]; ++e)
                for(f=0; f<lengths[p[5]-1]; ++f){
                  check=1;
                  for(i=0; i<clen; ++i){
                    sum=corr[p[0]-1][a][i]+
                      corr[p[1]-1][b][i]+
                      corr[p[2]-1][c][i]+
                      corr[p[3]-1][d][i]+
                      corr[p[4]-1][e][i]+
                      corr[p[5]-1][f][i];
                    if(sum%3){
                      --check;
                      break;
                    }
                  }
                  if(check){
                    G=permutations(z,zz,sub[p[0]-1][a],sub[p[1]-1][b],sub[p[2]-1][c],sub[p[3]-1][d],sub[p[4]-1][e],sub[p[5]-1][f]);
                    if(G!=NULL && D==min_wt(G)){
                      for(i=len; i<Len; ++i)
                        printf("%2lu",G[0][i]);
                      printf("\n");
                    }
                    if(G!=NULL){
                      for(i=0; i<len; ++i)
                        free(G[i]);
                      free(G);
                    }
                  }
                }
    }
    if(k==6 && p[0]<=LEN){
      for(a=0; a<lengths[p[0]-1]; ++a)
        for(b=0; b<lengths[p[1]-1]; ++b)
          for(c=0; c<lengths[p[2]-1]; ++c)
            for(d=0; d<lengths[p[3]-1]; ++d)
              for(e=0; e<lengths[p[4]-1]; ++e)
                for(f=0; f<lengths[p[5]-1]; ++f)
                  for(g=0; g<lengths[p[6]-1]; ++g){
                    check=1;
                    for(i=0; i<clen; ++i){
                      sum=corr[p[0]-1][a][i]+
                        corr[p[1]-1][b][i]+
                        corr[p[2]-1][c][i]+
                        corr[p[3]-1][d][i]+
                        corr[p[4]-1][e][i]+
                        corr[p[5]-1][f][i]+
                        corr[p[6]-1][g][i];
                      if(sum%3){
                        --check;
                        break;
                      }
                    }
                    if(check){
                      G=permutations(z,sub[p[0]-1][a],sub[p[1]-1][b],sub[p[2]-1][c],sub[p[3]-1][d],sub[p[4]-1][e],sub[p[5]-1][f],sub[p[6]-1][g]);
                      if(G!=NULL && D==min_wt(G)){
                        for(i=len; i<Len; ++i)
                          printf("%2lu",G[0][i]);
                        printf("\n");
                      }
                      if(G!=NULL){
                        for(i=0; i<len; ++i)
                          free(G[i]);
                        free(G);
                      }
                    }
                  }
    }
    if(k==7 && p[0]<=LEN){
      for(a=0; a<lengths[p[0]-1]; ++a)
        for(b=0; b<lengths[p[1]-1]; ++b)
          for(c=0; c<lengths[p[2]-1]; ++c)
            for(d=0; d<lengths[p[3]-1]; ++d)
              for(e=0; e<lengths[p[4]-1]; ++e)
                for(f=0; f<lengths[p[5]-1]; ++f)
                  for(g=0; g<lengths[p[6]-1]; ++g)
                    for(h=0; h<lengths[p[7]-1]; ++h){
                      check=1;
                      for(i=0; i<clen; ++i){
                        sum=corr[p[0]-1][a][i]+
                          corr[p[1]-1][b][i]+
                          corr[p[2]-1][c][i]+
                          corr[p[3]-1][d][i]+
                          corr[p[4]-1][e][i]+
                          corr[p[5]-1][f][i]+
                          corr[p[6]-1][g][i]+
                          corr[p[7]-1][h][i];
                        if(sum%3){
                          --check;
                          break;
                        }
                      }
                      if(check){
                        G=permutations(sub[p[0]-1][a],sub[p[1]-1][b],sub[p[2]-1][c],sub[p[3]-1][d],sub[p[4]-1][e],sub[p[5]-1][f],sub[p[6]-1][g],sub[p[7]-1][h]);
                        if(G!=NULL && D==min_wt(G)){
                          for(i=len; i<Len; ++i)
                            printf("%2lu",G[0][i]);
                          printf("\n");
                        }
                        if(G!=NULL){
                          for(i=0; i<len; ++i)
                            free(G[i]);
                          free(G);
                        }
                      }
                    }
    }

    int rem_val=0;
    while(k>=0 && p[k]==1){
      ++rem_val;
      --k;
    }

    if(k<0) return;

    --p[k];
    ++rem_val;

    while(rem_val>p[k]){
      p[k+1]=p[k];
      rem_val-=p[k++];
    }
    p[++k]=rem_val;
  }
  free(z);
  free(zz);
  free(zzz);
}

// driver ////////////////////////////////////////////////////////////////////////

int main(void)
{
  size_t i,j,len=8*LEN,Len=16*LEN;

  UINT lengths[LEN];
  {
    UINT c[LEN+1][LEN+1];
    for(i=0; i<=LEN; ++i)
      for(j=0; j<=i; ++j){
        if(j==0 || j==i) c[i][j]=(UINT)1;
        else c[i][j]=c[i-1][j]+c[i-1][j-1];
      }
    for(i=1; i<=LEN; lengths[i-1]=(((UINT)1)<<i)*c[LEN][i],++i);
  }

# if LEN%2
  unsigned clen=(LEN-1)/2;
# else

# if NEGA
  unsigned clen=(LEN-2)/2;
# else
  unsigned clen=LEN/2;
# endif

# endif

  UINT **sub[LEN],
    **corr[LEN];
  {
    size_t index;
    UINT counter,lim,num;
    unsigned check;
    for(size_t k=1; k<=LEN; ++k){
      sub[k-1]=(UINT**)malloc(lengths[k-1]*sizeof(UINT*));
      corr[k-1]=(UINT**)malloc(lengths[k-1]*sizeof(UINT*));
      for(i=0; i<lengths[k-1]; ++i){
        sub[k-1][i]=(UINT*)calloc(LEN,sizeof(UINT));
        corr[k-1][i]=(UINT*)malloc(clen*sizeof(UINT));
      }
      index=0;
      lim=(((UINT)1)<<k);
      UINT dat[k];
      for(i=0; i<k; dat[i]=(UINT)i,++i);
      for(counter=0; counter<lim; ++counter,++index){
        for(i=0; i<k; ++i){
          if(counter & (((UINT)1)<<i)) sub[k-1][index][dat[i]]=2;
          else sub[k-1][index][dat[i]]=1;
        }
        corrs(sub[k-1][index],corr[k-1][index],clen);
      }

      for(check=1; check;){
        for(i=0; i<k; ++i){
          num=dat[i]+1;
          if(!binsearch(dat,i+1,k-1,num)){
            if(i==k-1 && num==LEN){
              --check;
              break;
            }
            ++dat[i];
            for(j=0; j<i; dat[j]=(UINT)j,++j);
            for(counter=0; counter<lim; ++counter,++index){
              for(j=0; j<k; ++j){
                if(counter & (((UINT)1)<<j)) 
                  sub[k-1][index][dat[j]]=2;
                else
                  sub[k-1][index][dat[j]]=1;
              }
              corrs(sub[k-1][index],corr[k-1][index],clen);
            }
            break;
          }
        }
      }
    }
  }

# if LEN%3==2
  {
    size_t a,b,c,d;
    UINT *z=(UINT*)calloc(LEN,sizeof(UINT)),
      *zz=(UINT*)calloc(LEN,sizeof(UINT)),
      *zzz=(UINT*)calloc(LEN,sizeof(UINT)),
      *zzzz=(UINT*)calloc(LEN,sizeof(UINT)),
      **G,sum;
    unsigned check,D=3*floor(4*LEN/3.0)+3;
    for(a=0; a<lengths[LEN-1]; ++a)
      for(b=0; b<lengths[LEN-1]; ++b)
        for(c=0; c<lengths[LEN-1]; ++c)
          for(d=0; d<lengths[LEN-1]; ++d){
            check=1;
            for(i=0; i<clen; ++i){
              sum=corr[LEN-1][a][i]+
                corr[LEN-1][b][i]+
                corr[LEN-1][c][i]+
                corr[LEN-1][d][i];
              if(sum%3){
                --check;
                break;
              }
            }
            if(check){
              G=permutations(z,zz,zzz,zzzz,sub[LEN-1][a],sub[LEN-1][b],sub[LEN-1][c],sub[LEN-1][d]);
              if(G!=NULL && D==min_wt(G)){
                for(i=len; i<Len; ++i)
                  printf("%2lu",G[0][i]);
                printf("\n");
              }
              if(G!=NULL){
                for(i=0; i<len; ++i)
                  free(G[i]);
                free(G);
              }
            }
          }
    free(z);
    free(zz);
    free(zzz);
    free(zzzz);
  }
# endif

  int d=3*floor(2*LEN/3.0)+2;
  for(int w=d; w<=len; ++w)
    if(w%3 == 2)
      partitions(w,sub,corr,lengths,clen);

  for(i=0; i<LEN; ++i){
    for(j=0; j<lengths[i]; ++j){
      free(sub[i][j]);
      free(corr[i][j]);
    }
    free(sub[i]);
    free(corr[i]);
  }

  exit(EXIT_SUCCESS);
}
