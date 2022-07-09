# include "./search2.h"

typedef uint_fast32_t UINT;
typedef uint64_t MASK;

// misc ////////////////////////////////////////////////////////////////////////

void **transpose(size_t nrows, size_t ncols,
                 UINT *restrict A[restrict static 1],
                 UINT *restrict B[restrict static 1])
{
  size_t i,j;
  for(i=0; i<ncols; ++i)
    for(j=0; j<nrows; B[i][j]=A[j][i],++j);
}

static inline
void copy_mat(size_t len, size_t Len,
              UINT *restrict A[restrict static 1],
              UINT *restrict B[restrict static 1])
{
  size_t i,j;
  for(i=0; i<Len; ++i)
    for(j=0; j<len; ++j)
      B[i][j]=A[i][j];
}

// generator matrix /////////////////////////////////////////////////////////////

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

UINT **tcirc(UINT a[restrict static 1], unsigned minus)
{
    size_t i,j;
    UINT **A=(UINT**)malloc(LEN*sizeof(UINT*));
    for(i=0; i<LEN; A[i++]=(UINT*)malloc(LEN*sizeof(UINT)));
    if(minus){
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
    else{
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
    return A;
}

UINT **gen_mat(UINT a[restrict static 1],
               UINT b[restrict static 1])
{
    size_t i,j,len=2*LEN,Len=4*LEN,q,r;
    UINT **A=circ(a),
         **AT=tcirc(a,0),
         **B=circ(b),
         **BT=tcirc(b,1),
         **P[DIM][DIM]={{A,B},{BT,AT},},
         **G=(UINT**)malloc(len*sizeof(UINT*));
    for(i=0; i<len; G[i++]=(UINT*)malloc(Len*sizeof(UINT)));

    for(i=0; i<len; ++i)
        for(j=0; j<len; ++j){
            if(i == j) G[i][j]=1;
            else G[i][j]=0;
        }
    for(i=0; i<len; ++i){
        q=(i-(i%LEN))/LEN;
        r=i%LEN;
        for(j=0; j<len; ++j)
            G[i][len+j]=P[q][(j-(j%LEN))/LEN][r][j%LEN];
    }
    for(i=0; i<LEN; ++i){
        free(A[i]); free(AT[i]);
        free(B[i]); free(BT[i]);
    }
    free(A); free(AT);
    free(B); free(BT);

    return G;
}

// partitions ///////////////////////////////////////////////////////////////////

typedef struct {
  size_t index;
  MASK A[2*LEN];
} partition_t;

static inline
void partition_init(size_t len, partition_t p[static 1])
{
  p->index=0;
  for(size_t i=0; i<len; p->A[i]=0, ++i);
}

static inline
void partition_free(partition_t p[static 1])
{
  p->index=0;
}

static inline
void set_union(MASK *a, MASK *b)
{
  *a|=*b;
  *b=0;
}

static inline
void partition_reset(partition_t p[static 1])
{
  for(size_t i=0; i<p->index; ++i)
    if(!(p->A[i])){
      p->A[i]=p->A[--(p->index)];
      --i;
    }
}

static inline
size_t partition_index(size_t j, partition_t p[static 1])
{
  for(size_t i=0; i<p->index; ++i)
    if(p->A[i] & (1ULL << j))
      return i;
  return SIZE_MAX;
}

static inline
unsigned partition_contains(size_t j, MASK p)
{
  if(p & (1ULL << j)) return 1U;
  return 0U;
}

// matrices /////////////////////////////////////////////////////////////////////

static inline
UINT **matrix_init(size_t len, size_t Len)
{
  UINT **A=(UINT**)malloc(Len*sizeof(UINT*));
  for(size_t i=0; i<Len; A[i++]=(UINT*)malloc(len*sizeof(UINT)));
  return A;
}

static inline
void matrix_zeros(size_t len, size_t Len, UINT *A[static 1])
{
  size_t i,j;
  for(i=0; i<Len; ++i)
    for(j=0; j<len; ++j)
      A[i][j]=0;
}

void matrix_free(size_t Len, UINT *A[static 1])
{
  if(A == NULL) return;
  for(size_t i=0; i<Len; free(A[i]),++i);
  free(A);
  A=NULL;
}

unsigned matrix_equal(size_t len, size_t Len,
                      UINT *restrict A[restrict static 1],
                      UINT *restrict B[restrict static 1])
{
  size_t i,j;
  for(i=0; i<Len; ++i)
    for(j=0; j<len; ++j)
      if(A[i][j] != B[i][j])
        return 0U;
  return 1U;

}

static inline
void matrix_col_mult(size_t j, size_t Len, UINT *A[static 1])
{
  for(size_t i=0; i<Len; A[i][j]=(2*A[i][j])%3,++i);
}

static inline
void matrix_row_mult(size_t j, size_t len, UINT *A[static 1])
{
  for(size_t i=0; i<len; A[j][i]=(2*A[j][i])%3,++i);
}

MASK matrix_row_support(size_t j, size_t len, UINT *A[static 1])
{
  MASK row=0;
  for(size_t i=0; i<len; ++i)
    if(A[j][i])
      row |= 1ULL<<i;
  return row;
}

UINT **generator(UINT a[restrict static 1],
                 UINT b[restrict static 1])
{
  size_t i,len=2*LEN,Len=4*LEN;
  UINT **A=gen_mat(a,b),
    **B=matrix_init(len,Len);
  transpose(len,Len,A,B);
  matrix_free(len,A);
  return B;
}

// gaussian eliminiation of col /////////////////////////////////////////////////

static inline
void swap_rows(UINT **A, size_t r1, size_t r2)
{
  UINT *temp=A[r1];
  A[r1]=A[r2];
  A[r2]=temp;
}

static inline
void row_mult(UINT *A[static 1], size_t ncols, size_t r, UINT mult)
{
  for(size_t i=0; i<ncols; A[r][i]=(mult*A[r][i])%3,++i);
}

static inline
void row_mult_add(UINT *A[static 1], size_t ncols, size_t r1, size_t r2, UINT mult)
{
  for(size_t i=0; i<ncols; ++i)
    A[r2][i]=(A[r2][i] + (mult*A[r1][i]))%3;
}

void gauss_elim(size_t k, size_t len, size_t Len, size_t rank,
                UINT *restrict m[restrict static 1])
{
  size_t i;
  UINT **A=matrix_init(Len,len),mult;

  transpose(Len,len,m,A);

  if(!A[rank][k]){
    for(i=rank+1; i<len; ++i)
      if(A[i][k]) break;
    swap_rows((UINT**)A,rank,i);
  }

  if(A[rank][k] == 2) row_mult(A,Len,rank,2);

  for(i=rank+1; i<len; ++i)
    if(A[i][k]){
      mult=3-A[i][k];
      row_mult_add((UINT**)A,Len,rank,i,mult);
    }

  for(i=0; i<rank; ++i)
    if(A[i][k]){
      mult=3-A[i][k];
      row_mult_add((UINT**)A,Len,rank,i,mult);
    }

  transpose(len,Len,A,m);
  matrix_free(len,(UINT**)A);
}

// order functions ///////////////////////////////////////////////////////////

int colex(size_t len,
          UINT a[restrict static 1],
          UINT b[restrict static 1])
{
  for(size_t i=len; i; --i){
    if(a[i-1] == b[i-1]) continue;
    if(a[i-1] < b[i-1]) return -1;
    return 1;
  }
  return 0;
}

int lex(size_t i,
        size_t len,
        UINT *restrict A[restrict static 1],
        UINT *restrict B[restrict static 1])
{
  int check;
  for(size_t j=0; j<i; ++j){
    if((check=colex(len,A[j],B[j])) == 0)
      continue;
    return check;
  }
  return 0;
}

// support signature ///////////////////////////////////////////////////////////

MASK *support_generate(size_t i, size_t len, size_t Len, UINT **A)
{
  size_t j;
  MASK span=0,row,*supports=(MASK*)malloc(Len*sizeof(MASK));

  for(j=0; j<i; ++j){
    supports[j]=matrix_row_support(j,len,A);
    span |= supports[j];
  }

  for(size_t j=i; j<Len; ++j){
    row=matrix_row_support(j,len,A);
    if(row == (row & span))
      supports[j]=row;
    else
      supports[j]=0ULL;
  }

  return supports;
}

int compare_support(void const *aa, void const *bb)
{
  MASK a=*(MASK*)aa, b=*(MASK*)bb;

  if(a == 0) return 1;
  if(b == 0) return -1;

  if(a < b) return -1;
  if(b < a) return 1;

  return 0;
}

unsigned *support_sort(size_t Len,
                       size_t nreps,
                       unsigned reps[static 1],
                       MASK supports[restrict static 1],
                       MASK candidate[restrict static 1])
{
  size_t i,j;
  unsigned *perm=(unsigned*)malloc(Len*sizeof(unsigned));

  for(i=0; i<Len; candidate[i]=supports[i],++i);

  for(i=0; i<nreps-1; ++i)
    qsort(candidate+i,reps[i+1]-reps[i],sizeof(MASK),compare_support);

  for(i=0; i<Len; ++i)
    for(j=0; j<Len; ++j)
      if(supports[i] == candidate[j]){
        perm[i]=j;
        break;
      }

  return perm;
}

unsigned *support_partition(size_t Len,
                            size_t nreps,
                            unsigned reps[static 1],
                            size_t new_nreps[static 1],
                            MASK candidate[static 1])
{
  size_t i,j,m,index;
  unsigned *new_reps=(unsigned*)calloc((Len+1),sizeof(unsigned));

  for(i=0,index=0; i<nreps-1; ++i){
    new_reps[index++]=reps[i];
    m=reps[i+1]-reps[i];
    for(j=1; j<m; ++j){
      if(candidate[j-1] != candidate[j]){
        new_reps[index++]=reps[i]+j;
        continue;
      }
    }
  }
  new_reps[index++]=Len;
  *new_nreps=index;

  return new_reps;
}

// canonization ////////////////////////////////////////////////////////////////

void RC(size_t lam, size_t len, size_t Len,
        UINT *m[static 1], partition_t P[static 1])
{
  size_t l;
  MASK p=P->A[lam],row;

  for(l=0; l<len; ++l)
    if(partition_contains(l,p))
      matrix_col_mult(l,Len,m);

  for(l=0; l<Len; ++l){
    row=matrix_row_support(l,len,m);
    if(row == (row & p))
      matrix_row_mult(l,len,m);
  }
}

void semicanonical(size_t i, size_t len, size_t Len,
                   UINT *restrict m[restrict static 1],
                   size_t rank[static 1],
                   partition_t P[static 1],
                   UINT *restrict ith_m[restrict static 1])
{
  if(i == 0){
    copy_mat(len,Len,m,ith_m);
    /* *rank=0; */
    /* partition_init(len,P); */
    return;
  }

  semicanonical(i-1,len,Len,m,rank,P,ith_m);

  size_t j,union_index,lam;
  unsigned check=0;

  for(j=*rank; j<len; ++j)
    if(m[i-1][j]){
      check=1;
      break;
    }

  if(check){
    gauss_elim(i-1,len,Len,*rank,ith_m);
    P->A[P->index++]=1ULL<<*rank;
    (*rank)++;
  }
  else{
    union_index=SIZE_MAX;

    for(j=*rank; j; --j){
      if(ith_m[i-1][j-1] == 0)
        continue;

      lam=partition_index(j-1,P);
      assert(lam != SIZE_MAX);

      if(union_index == SIZE_MAX){
        union_index=lam;
        if(ith_m[i-1][j-1] == 2)
          matrix_row_mult(i-1,len,(UINT**)ith_m);
      }
      else if(lam != union_index){
        if(ith_m[i-1][j-1] == 2)
          RC(lam,len,Len,(UINT**)ith_m,P);
        set_union(&P->A[union_index],&P->A[lam]);
      }
    }
  }
  partition_reset(P);
}

void canonical_inner(size_t i, size_t len, size_t Len,
                     unsigned cand_check[static 1],
                     UINT *restrict candidate[restrict static 1],
                     UINT *restrict m[restrict static 1])
{
  int check;
  size_t rank;
  partition_t P;
  UINT **gamma=matrix_init(len,Len);
  matrix_zeros(len,Len,gamma);

  for(size_t j=i-1; j<Len; ++j){
    rank=0;
    partition_init(len,&P);

    swap_rows((UINT**)m,i-1,j);
    semicanonical(i,len,Len,m,&rank,&P,gamma);
    swap_rows((UINT**)m,i-1,j);

    check=lex(i,len,gamma,candidate);

    if(*cand_check == 1 && check == 1)
      continue;
    if(*cand_check == 1 && check == -1)
      *cand_check=0;
    if(i<Len)
      canonical_inner(i+1,len,Len,cand_check,candidate,gamma);
    else if(*cand_check == 0){
      copy_mat(len,Len,gamma,candidate);
      *cand_check=1;
    }
  }

  matrix_free(Len,gamma);
}

UINT **canonical(size_t len, size_t Len,
                 UINT *restrict m[restrict static 1])
{
  unsigned cand_check=0;
  UINT **candidate=matrix_init(len,Len);
  matrix_zeros(len,Len,candidate);
  canonical_inner(1,len,Len,&cand_check,candidate,m);
  return candidate;
}

// driver ///////////////////////////////////////////////////////////////////////

int main(void)
{
  size_t i,j,len=2*LEN,Len=4*LEN;
  UINT a[]={1,1,},
    b[]={1,1,},
    **G=generator(a,b),**H;

  H=canonical(len,Len,G);

  for(i=0; i<len; ++i){
    for(j=0; j<Len; printf("%2lu",G[j++][i]));
    printf("\n");
  }
  printf("\n");

  MASK *sup=support_generate(len,len,Len,G);
  sup[0]+=2;

  for(i=0; i<len; ++i){
    for(j=0; j<Len; ++j){
      if(sup[j] & (1ULL << i)) printf("1 ");
      else printf("0 ");
    }
    printf("\n");
  }
  printf("\n");

  unsigned reps[]={0,Len};
  size_t nreps=2;
  MASK *cand=(MASK*)malloc(Len*sizeof(MASK));
  unsigned *perm=support_sort(Len,nreps,reps,sup,cand);

  for(i=0; i<len; ++i){
    for(j=0; j<Len; ++j){
      if(cand[j] & (1ULL << i)) printf("1 ");
      else printf("0 ");
    }
    printf("\n");
  }
  printf("\n");

  for(i=0; i<Len; printf("%-2zu",i++));
  printf("\n");
  for(i=0; i<Len; printf("%-2u",perm[i++]));
  printf("\n\n");

  size_t new_nreps;
  unsigned *new_reps=support_partition(Len,nreps,reps,&new_nreps,cand);

  for(i=0; i<Len+1; printf("%-3u",new_reps[i++]));
  printf("\n");
  printf("%zu\n",new_nreps);

  free(new_reps);
  free(perm);
  free(cand);
  free(sup);
  matrix_free(Len,H);
  matrix_free(Len,G);

  exit(EXIT_SUCCESS);
}
