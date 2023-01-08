#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include <omp.h>

long numThreads=8;
long dim=1000;


struct Matrix{
    unsigned int x;
    unsigned int y;
    double *m; 
};


struct Matrix CreateMatrix(unsigned int x,unsigned int y){
    struct Matrix lin;
    lin.x=x;
    lin.y=y;
    lin.m=(double *)malloc(x*y*sizeof(double));
    int i;
    for(i=0;i<x*y;i++){
        lin.m[i]=0;
    }
    return lin;
}

void *PointMatrix(struct Matrix m){
    int i,j;
    for(i=0;i<m.x;i++){
        for(j=0;j<m.y;j++){
            printf("%f\t",m.m[i*m.y+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void*positive_definite_matrix(double*m){
    int i, j;
	double sum = 1;
	for (i = 0;i < dim;i++) {
		for (j = 0;j < dim;j++) {
			if (j < i)
				m[i * dim + j] = 0;
			else 
			    m[i * dim + j] = sum++;
		}
	}
}

void transpose(struct Matrix M,struct Matrix MT){
    int i,j;
    for(i=0;i<MT.x;i++){
        for(j=0;j<MT.y;j++){
            MT.m[i*MT.y+j]=M.m[j*M.y+i];
        }
    }
}

double vectorMul(struct Matrix mul1,struct Matrix mul2){
    int i;
    double sum=0;
    for(i=0;i<mul1.y;i++){
        sum+=mul1.m[i]*mul2.m[i];
    }
    return sum;
}

void mulReal(double r,struct Matrix A,struct Matrix rA){
    int i,j;
    for(i=0;i<rA.x;i++){
        for(j=0;j<rA.y;j++){
            rA.m[i*rA.y+j]=r*A.m[i*A.y+j];
        }
    }
}

void addMatrix(struct Matrix a1,struct Matrix a2,struct Matrix sum){
    int i,j;
    for(i=0;i<sum.x;i++){
        for(j=0;j<sum.y;j++){
            sum.m[i*sum.y+j]=a1.m[i*sum.y+j]+a2.m[i*sum.y+j];
        }
    }
}

void subMatrix(struct Matrix a1,struct Matrix a2,struct Matrix sum){
    int i,j;
    for(i=0;i<sum.x;i++){
        for(j=0;j<sum.y;j++){
            sum.m[i*sum.y+j]=a1.m[i*sum.y+j]-a2.m[i*sum.y+j];
        }
    }
}

void fuzhi(struct Matrix a1,struct Matrix a2){
    int i,j;
    for(i=0;i<a1.x;i++){
        for(j=0;j<a1.y;j++){
            a2.m[i*a1.y+j]=a1.m[i*a1.y+j];
        }
    }
}

double finish(struct Matrix r){
    int i,j;
    double sum=0;
    double lin;
    for(i=0;i<r.x;i++){
        for(j=0;j<r.y;j++){
            lin=fabs(r.m[i*r.y+j]);
            sum+=lin*lin;
        }
    }
    return sum;
}

void resetmatrix(struct Matrix M){
    long i,j;
    for(i=0;i<M.x;i++){
        for(j=0;j<M.y;j++){
            M.m[i*M.y+j]=0;
        }
    }

}
void MultMatrix(struct Matrix mul1,struct Matrix mul2,struct Matrix prod){
	long i,j,k;
#pragma omp parallel
{
#pragma omp for schedule(dynamic)
    for(i=0;i<prod.x;i++){
		for(j=0;j<prod.y;j++){
			for(k=0;k<mul1.y;k++){
				prod.m[i*prod.y+j]+=mul1.m[i*mul1.y+k]*mul2.m[k*mul2.y+j];
			}
		}
	}
}
#pragma omp barrier
}


int main()
{

	int i,j,k;
    struct timeval t1,t2;

	struct Matrix A=CreateMatrix(dim,dim);
    positive_definite_matrix(A.m);
    struct Matrix X=CreateMatrix(dim,1);
    for(i=0;i<dim;i++){
        X.m[i]=0;
    }
    struct Matrix B=CreateMatrix(dim,1);
    for(i=0;i<dim;i++){
        B.m[i]=2;
    }
	gettimeofday(&t1,NULL);
    omp_set_num_threads(numThreads);
	struct Matrix Ax,r1,rt,p1,pt,pa,ap1,x2,aA,aAp,r2,r2t,bp,p2;
    double al,r1r,pap,fr,r2r,bat;
	Ax=CreateMatrix(dim,1);
    r1=CreateMatrix(dim,1);
    rt=CreateMatrix(1,dim);
    p1=CreateMatrix(dim,1);
    pt=CreateMatrix(1,dim);
    pa=CreateMatrix(1,dim);
    ap1=CreateMatrix(dim,1);
    x2=CreateMatrix(dim,1);
    aA=CreateMatrix(dim,dim);
    aAp=CreateMatrix(dim,1);
    r2=CreateMatrix(dim,1);
    r2t=CreateMatrix(dim,1);
    bp=CreateMatrix(dim,1);
    p2=CreateMatrix(dim,1);

    int ii,jj,kk;
    int sum=0;

//     double a1[dim][dim];
//     double a2[dim][dim];
//     double a3[dim][dim];
//     for(i=0;i<dim;i++){
//         for(j=0;j<dim;j++){
//             a1[i][j]=1;
//             a2[i][j]=2;
//             a3[i][j]=0;
//         }
//     }
// #pragma omp parallel shared(a1,a2,a3) private(ii,jj,kk)
// {
//    #pragma omp for schedule(dynamic)
//     for(ii=0;ii<Ax.x;ii++){
// 		for(jj=0;jj<Ax.y;jj++){
// 		    for(kk=0;kk<A.y;kk++){
// 			    a3[ii][jj]+=a2[ii][kk]*a3[kk][jj];
// 		    }
// 		}
// 	}
// }
//#pragma omp barrier

    MultMatrix(A,X,Ax);
	subMatrix(B,Ax,r1);
    fuzhi(r1,p1);
    i=0;
	while(i<20){
        transpose(r1,rt);

        r1r=vectorMul(rt,r1);
        transpose(p1,pt);

        resetmatrix(pa);
        MultMatrix(pt,A,pa);

        // #pragma omp parallel for
        // for(ii=0;ii<pa.x;ii++){
		//     for(jj=0;jj<pa.y;jj++){
		// 	    for(kk=0;kk<pt.y;kk++){
		// 		    pa.m[ii*pa.y+jj]+=pt.m[ii*pt.y+kk]*A.m[kk*A.y+jj];
		// 	    }
		//     }
	    // }
        // #pragma omp barrier
        pap=vectorMul(pa,p1);
 
        al=r1r/pap;

        mulReal(al,p1,ap1);

        addMatrix(X,ap1,x2);

        mulReal(al,A,aA);

        resetmatrix(aAp);
        MultMatrix(aA,p1,aAp);
        // #pragma omp parallel for
        // for(ii=0;ii<aAp.x;ii++){
		//     for(jj=0;jj<aAp.y;jj++){
		// 	    for(kk=0;kk<aA.y;kk++){
		// 		    aAp.m[ii*aAp.y+jj]+=aA.m[ii*aA.y+kk]*p1.m[kk*p1.y+jj];
		// 	    }
		//     }
	    // }
        // #pragma omp barrier
        subMatrix(r1,aAp,r2);

        fr=finish(r2);


        if(fr<0.001)
            break;
        
        transpose(r2,r2t);
        r2r=vectorMul(r2t,r2);
        bat=r2r/r1r;
        mulReal(bat,p1,bp);
        addMatrix(r2,bp,p2);

        fuzhi(p2,p1);
        fuzhi(r2,r1);
        fuzhi(x2,X);
        i++;
    }
    
    gettimeofday(&t2,NULL);
    printf("Количество итерации:%d\n",i);
    double time=t2.tv_sec-t1.tv_sec+(double)(t2.tv_usec-t1.tv_usec)/1000000;
    printf("Размерность матрицы:%d\n",dim);
    printf("Количество потоков:%d\nВремя:%fs\n",numThreads,time);
    
    return 0;
}