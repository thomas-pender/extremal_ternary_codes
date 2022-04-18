# include "./search4.h"

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
unsigned min(unsigned a, unsigned b)
{
  if(a <= b) return a;
  else return b;
}

static inline
void copy_mat(UINT *restrict A[restrict static 1],
              UINT *restrict B[restrict static 1])
{
  size_t i,j,len=4*LEN,Len=8*LEN;
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

void transpose(UINT *restrict A[restrict static 1],
               UINT *restrict B[restrict static 1])
{
  size_t i,j,len=4*LEN,Len=8*LEN;
  for(i=0; i<Len; ++i)
    for(j=0; j<len; B[i][j]=A[j][i],++j);
}

void corrs(UINT a[restrict static 1], 
           UINT corr[restrict static 1], 
           unsigned len)
{
  size_t i,j;
  UINT sum;
  for(j=1; j<=len; ++j){
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

// generator matrix //////////////////////////////////////////////////////////////

UINT **circ(UINT a[restrict static 1])
{
  size_t i,j;
  UINT **A=(UINT**)malloc(LEN*sizeof(UINT*));
  for(i=0; i<LEN; A[i++]=(UINT*)malloc(LEN*sizeof(UINT)));
  for(i=0; i<LEN; A[0][i]=a[i],++i);
  for(i=1; i<LEN; ++i)
    for(j=0; j<LEN; ++j){
      if(j) A[i][j]=A[i-1][j-1];
      else{
# if NEGA
        A[i][j]=(2*A[i-1][LEN-1])%3;
# else
        A[i][j]=A[i-1][LEN-1];
# endif
      }
    }
  return A;
}

UINT **back_circ(UINT a[restrict static 1], unsigned transpose)
{
  size_t i,j;
  UINT **A=(UINT**)malloc(LEN*sizeof(UINT*));
  for(i=0; i<LEN; A[i++]=(UINT*)malloc(LEN*sizeof(UINT)));
  if(transpose){
    A[0][LEN-1]=a[0];
# if NEGA
    for(i=0; i<LEN-1; A[0][i]=(2*a[i+1])%3,++i);
# else
    for(i=0; i<LEN-1; A[0][i]=a[i+1],++i);
# endif
  }
  else
    for(i=0; i<LEN; A[0][i]=a[LEN-i-1],++i);
  for(i=1; i<LEN; ++i)
    for(j=0; j<LEN; ++j){
      if(j == LEN-1){
# if NEGA
        A[i][j]=(2*A[i-1][0])%3;
# else
        A[i][j]=A[i-1][0];
# endif
      }
      else
        A[i][j]=A[i-1][j+1];
    }
  return A;
}

UINT **minus(UINT *restrict A[restrict static 1])
{
  size_t i,j;
  UINT **B=(UINT**)malloc(LEN*sizeof(UINT*));
  for(i=0; i<LEN; ++i){
    B[i]=(UINT*)malloc(LEN*sizeof(UINT));
    for(j=0; j<LEN; B[i][j]=(2*A[i][j])%3,++j);
  }
  return B;
}

UINT **gen_mat(UINT a[restrict static 1],
               UINT b[restrict static 1],
               UINT c[restrict static 1],
               UINT d[restrict static 1])
{
  size_t i,j,len=4*LEN,Len=8*LEN;

  UINT **A=circ(a),
    **B=back_circ(b,0),
    **MB=minus(B),
    **BT=back_circ(b,1),
    **MBT=minus(BT),
    **C=back_circ(c,0),
    **MC=minus(C),
    **CT=back_circ(c,1),
    **MCT=minus(CT),
    **D=back_circ(d,0),
    **MD=minus(D),
    **DT=back_circ(d,1),
    **MDT=minus(DT),
    **P[DIM][DIM]={{A,B,C,D,},{MB,A,MDT,CT,},{MC,DT,A,MBT,},{MD,MCT,BT,A,},};
  UINT **G=(UINT**)malloc(len*sizeof(UINT*));
  for(i=0; i<len; G[i++]=(UINT*)malloc(Len*sizeof(UINT)));
  for(i=0; i<len; ++i)
    for(j=0; j<len; ++j){
      if(i == j) G[i][j]=1;
      else G[i][j]=0;
    }
  for(i=0; i<len; ++i)
    for(j=0; j<len; ++j)
      G[i][len+j]=P[(i-(i%LEN))/LEN][(j-(j%LEN))/LEN][i%LEN][j%LEN];

  for(i=0; i<LEN; ++i){
    free(A[i]);
    free(B[i]); free(MB[i]); free(BT[i]); free(MBT[i]);
    free(C[i]); free(MC[i]); free(CT[i]); free(MCT[i]);
    free(D[i]); free(MD[i]); free(DT[i]); free(MDT[i]);
  }
  free(A);
  free(B); free(MB); free(BT); free(MBT);
  free(C); free(MC); free(CT); free(MCT);
  free(D); free(MD); free(DT); free(MDT);

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
  unsigned check=1,len=4*LEN;
  for(size_t i=0; i<ncols; ++i){
    A[r2][i]=(A[r2][i] + (mult*A[r1][i]))%3;
    if(i>=len && check && A[r2][i]) --check;
  }
  return check;
}

UINT **rref(UINT *restrict G[restrict static 1])
{
  size_t i,j,k,nrows=4*LEN,ncols=8*LEN,row,col;
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
  size_t i,j,len=4*LEN,Len=8*LEN;
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
  size_t i,j,len=4*LEN;
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
  size_t i,len=4*LEN,Len=8*LEN;
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

// partitions ////////////////////////////////////////////////////////////////////

void partitions(int n,
                UINT *restrict*restrict sub[restrict static 1],
                UINT *restrict*restrict corr[restrict static 1],
                UINT lengths[restrict static 1],
                unsigned clen)
{
  register size_t i,a,b,c,d,len=4*LEN;
  unsigned check,D=3*floor(2*LEN/3.0)+3;
  UINT sum,**G,*z=(UINT*)calloc(LEN,sizeof(UINT));

  int p[n];
  int k=0;
  p[0]=n;

  while(1){
    if(k==2 && p[0]<=LEN){
      for(a=0; a<lengths[p[0]-1]; ++a)
        for(b=0; b<lengths[p[1]-1]; ++b)
          for(c=0; c<lengths[p[2]-1]; ++c){
            check=1;
            for(i=0; i<clen; ++i){
              sum=corr[p[0]-1][a][i]+
                corr[p[1]-1][b][i]+
                corr[p[2]-1][c][i];
              if(sum%3){
                --check;
                break;
              }
            }
            if(check){
              G=gen_mat(z,
                        sub[p[0]-1][a],
                        sub[p[1]-1][b],
                        sub[p[2]-1][c]);
              if(D == min_wt(G)){
                for(i=len; i<2*len; printf("%2lu",G[0][i++]));
                printf("\n");
                fflush(stdout);
              }
              for(i=0; i<len; free(G[i]),++i);
              free(G);
            }
          }
    }
    if(k==3 && p[0]<=LEN){
      for(a=0; a<lengths[p[0]-1]; ++a)
        for(b=0; b<lengths[p[1]-1]; ++b)
          for(c=0; c<lengths[p[2]-1]; ++c)
            for(d=0; d<lengths[p[3]-1]; ++d){
              check=1;
              for(i=0; i<clen; ++i){
                sum=corr[p[0]-1][a][i]+
                  corr[p[1]-1][b][i]+
                  corr[p[2]-1][c][i]+
                  corr[p[3]-1][d][i];
                if(sum%3){
                  --check;
                  break;
                }
              }
              if(check){
                G=gen_mat(sub[p[0]-1][a],
                          sub[p[1]-1][b],
                          sub[p[2]-1][c],
                          sub[p[3]-1][d]);
                if(D == min_wt(G)){
                  for(i=len; i<2*len; printf("%2lu",G[0][i++]));
                  printf("\n");
                  fflush(stdout);
                }
                for(i=0; i<len; free(G[i]),++i);
                free(G);
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
}

// driver ////////////////////////////////////////////////////////////////////////

int main(void)
{
  size_t i,j,len=4*LEN,Len=8*LEN;

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

# if LEN%3==1
  {
    size_t a,b;
    UINT *z=(UINT*)calloc(LEN,sizeof(UINT)),
      *zz=(UINT*)calloc(LEN,sizeof(UINT)),
      **G,sum;
    unsigned check,d=3*floor(2*LEN/3.0)+3;
    for(a=0; a<lengths[LEN-1]; ++a)
      for(b=0; b<lengths[LEN-1]; ++b){
        check=1;
        for(i=0; i<clen; ++i){
          sum=corr[LEN-1][a][i]+
            corr[LEN-1][b][i];
          if(sum%3){
            --check;
            break;
          }
        }
        if(check){
          G=gen_mat(z,zz,sub[LEN-1][a],sub[LEN-1][b]);
          if(d == min_wt(G)){
            for(i=len; i<Len; printf("%2lu",G[0][i++]));
            printf("\n");
            fflush(stdout);
          }
          for(i=0; i<len; free(G[i]),++i);
          free(G);
        }
      }
    free(z);
    free(zz);
  }
# endif

  int d=3*floor(2*LEN/3.0)+2;
  for(int w=d; w<len; ++w)
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
