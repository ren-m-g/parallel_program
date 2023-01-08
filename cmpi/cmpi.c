#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "mpi.h"

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
/// @brief 打印矩阵
/// @param m 矩阵
/// @return null
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


/// @brief 生成正定矩阵
/// @param 矩阵数组 
/// @return null
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


/// @brief 矩阵转置函数
/// @param M 原矩阵
/// @param MT 转置结果
void transpose(struct Matrix M,struct Matrix MT){
    int i,j;
    for(i=0;i<MT.x;i++){
        for(j=0;j<MT.y;j++){
            MT.m[i*MT.y+j]=M.m[j*M.y+i];
            //printf("MT[%d]=%f,M[%d}=%f\n",i*MT.y+j,MT.m[i*MT.y+j],j*M.y+i,M.m[j*M.y+i]);
        }
    }
}

/// @brief 向量乘法
/// @param mul1 乘数1
/// @param mul2 乘数2
/// @return 结果
double vectorMul(struct Matrix mul1,struct Matrix mul2){
    int i;
    double sum=0;
    for(i=0;i<mul1.y;i++){
        sum+=mul1.m[i]*mul2.m[i];
    }
    return sum;
}

/// @brief 矩阵与实数的乘法
/// @param r 实数
/// @param A 矩阵
/// @param rA 结果
void mulReal(double r,struct Matrix A,struct Matrix rA){
    int i,j;
    for(i=0;i<rA.x;i++){
        for(j=0;j<rA.y;j++){
            rA.m[i*rA.y+j]=r*A.m[i*A.y+j];
        }
    }
}

/// @brief 矩阵加法
/// @param a1 加数 
/// @param a2 被加数
/// @param sum 和
void addMatrix(struct Matrix a1,struct Matrix a2,struct Matrix sum){
    int i,j;
    for(i=0;i<sum.x;i++){
        for(j=0;j<sum.y;j++){
            sum.m[i*sum.y+j]=a1.m[i*sum.y+j]+a2.m[i*sum.y+j];
        }
    }
}

/// @brief 矩阵减法
/// @param a1 被减数
/// @param a2 减数
/// @param sum 差
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

void resetMatrix(struct Matrix A){
    int i,j;
    for(i=0;i<A.x;i++){
        for(j=0;j<A.y;j++){
            A.m[i*A.y+j]=0;
        }
    }
}

void MultMatrix(int id,int numprocs,struct Matrix mul1,struct Matrix mul2,struct Matrix prod){
    resetMatrix(prod);
    int i,j,k,X,Y;
    long num=dim/numprocs+1;
    long indexProd=id*num;
    for(i=0;i<num;i++){
        if(indexProd+i==prod.x*prod.y)
            break;
        X=(indexProd+i)/prod.y;
        Y=(indexProd+i)%prod.y;
        double sum=0;
        for(j=0;j<mul1.y;j++){
            sum+=mul1.m[X*mul1.y+j]*mul2.m[j*mul2.y+Y];
        }
        prod.m[indexProd+i]=sum;
    }

    //MPI_Barrier(MPI_COMM_WORLD);
    double *m_buffer=(double*)malloc(prod.x*prod.y*sizeof(double));
    if (id != 0){
        MPI_Send(prod.m, prod.x*prod.y, MPI_DOUBLE, 0, 99,
            MPI_COMM_WORLD);
    }
    else{
        for (k = 1; k < numprocs; k++) {
            MPI_Recv(m_buffer, prod.x*prod.y, MPI_DOUBLE, k, 99,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for(i=0;i<prod.x;i++){
                for(j=0;j<prod.y;j++){
                    prod.m[i*prod.y+j]+=m_buffer[i*prod.y+j];
                }
            }
        }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(prod.m, prod.x*prod.y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}


void main(int argc, char* argv[])
{
    int numprocs, myid, source;
    MPI_Status status;
    char message[100];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);



    int i,j,k;
    struct timeval t1,t2;
    long num=dim/numprocs+1;

    struct Matrix A=CreateMatrix(dim,dim);
    positive_definite_matrix(A.m);

    struct Matrix X=CreateMatrix(dim,1);
    for(i=0;i<dim;i++){
        X.m[i]=0;
    }

    struct Matrix B=CreateMatrix(dim,1);
    for(i=0;i<dim;i++){
        B.m[i]=100;
    }

    gettimeofday(&t1,NULL);

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

    MultMatrix(myid,numprocs,A,X,Ax);

    subMatrix(B,Ax,r1);
    fuzhi(r1,p1);
    i=0;

    while(i<20){
        transpose(r1,rt);
        r1r=vectorMul(rt,r1);
        transpose(p1,pt);

        //1
        MultMatrix(myid,numprocs,pt,A,pa);

        pap=vectorMul(pa,p1);
        al=r1r/pap;

        mulReal(al,p1,ap1);
        addMatrix(X,ap1,x2);

        mulReal(al,A,aA);

        //2
        MultMatrix(myid,numprocs,aA,p1,aAp);

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

    //PointMatrix(X);
    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&t2,NULL);
    if(myid==0){
        printf("Количество итерации:%d\n",i);
        double Xn=B.m[B.x*B.y-1]/A.m[A.x*A.y-1];
        double Xn_1=(B.m[B.x*B.y-2]-A.m[(A.x-1)*A.y-1]*Xn)/(A.m[(A.x-1)*A.y-2]);
        printf("Метод Гаусса:x(n):%f,x(n-1)=%f\n",Xn,Xn_1);
        printf("Метод сопряженных градиентов:x(n)=%f,x(n-1)=%f\n",x2.m[x2.x*x2.y-1],x2.m[x2.x*x2.y-2]);
        double time=t2.tv_sec-t1.tv_sec+(double)(t2.tv_usec-t1.tv_usec)/1000000;
        printf("Размерность матрицы:%d\n",dim);
        printf("Количество потоков:%d\nВремя:%fs\n",numprocs,time);


        //double time=t2.tv_sec-t1.tv_sec+(double)(t2.tv_usec-t1.tv_usec)/1000000;
        //printf("时间:%fs\n",time);
    }

    MPI_Finalize();
} /* end main */