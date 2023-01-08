#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
long numThreads=1;
long dim=1000;

void *PrintHello()
{
   printf("Hello World! It's me, thread #!\n");
}


struct Matrix{
    unsigned int x;
    unsigned int y;
    double *m; 
};

struct DataMultMatrix{
    unsigned int threadid;
    struct Matrix mul1;
    struct Matrix mul2;
    struct Matrix prod;
};

void initData(int sizedata,struct DataMultMatrix*data,struct Matrix mul1,struct Matrix mul2,struct Matrix prod){
    int i;
    for(i=0;i<sizedata;i++){
        data[i].mul1=mul1;
        data[i].mul2=mul2;
        data[i].prod=prod;
        data[i].threadid=0;
    }
}
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


void * MultMatrix(void * d) {
    struct DataMultMatrix data=*((struct DataMultMatrix*)d);
    unsigned int id = data.threadid;
    struct Matrix mul1 = data.mul1;
    struct Matrix mul2 = data.mul2;
    struct Matrix prod = data.prod;
    long num=(prod.x*prod.y)/numThreads+1;
    
    //printf("threadId=%d\n",id);

    long indexProd=num*id;
    int i,j,X,Y;
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

    return NULL;
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
void* my_thread(int numThreads,pthread_t* threads,struct DataMultMatrix*data){
    int rc;
    int i;
    for(i=0;i<numThreads;i++){
        data[i].threadid=i;
        rc = pthread_create(&threads[i], NULL, MultMatrix, (void *)&(data[i]));
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }  
    return NULL;
}
void* my_join(pthread_t* threads,int numThreads){
    int rc;
    int t;
    for(t=0;t<numThreads;t++)
    {
        rc = pthread_join(threads[t], NULL);
        if (rc != 0) {
            printf("1：等待线程失败");
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


int main(int argc, char *argv[])
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
        B.m[i]=10;
    }

    gettimeofday(&t1,NULL);
    pthread_t* threads=(pthread_t *)malloc(numThreads*sizeof(pthread_t));
    struct Matrix Ax,r1,rt,p1,pt,pa,ap1,x2,aA,aAp,r2,r2t,bp,p2;
    double al,r1r,pap,fr,r2r,bat;
    struct DataMultMatrix*data1=(struct DataMultMatrix*)malloc(numThreads*sizeof(struct DataMultMatrix));
    struct DataMultMatrix*data2=(struct DataMultMatrix*)malloc(numThreads*sizeof(struct DataMultMatrix));
    struct DataMultMatrix*data3=(struct DataMultMatrix*)malloc(numThreads*sizeof(struct DataMultMatrix));
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


    initData(numThreads,data1,A,X,Ax);
    my_thread(numThreads,threads,data1);
    my_join(threads,numThreads);
    subMatrix(B,Ax,r1);
    fuzhi(r1,p1);
    i=0;
    while(i<20){
        //printf("i=%d\n",i);
        transpose(r1,rt);
        //printf("r1\n");
        //PointMatrix(r1);
        //printf("rt\n");
        //PointMatrix(rt);


        r1r=vectorMul(rt,r1);
        //printf("r1r=%f\n",r1r);

        transpose(p1,pt);
        //printf("p1\n");
        //PointMatrix(p1);
        //printf("pt\n");
        //PointMatrix(pt);


        initData(numThreads,data2,pt,A,pa);
        my_thread(numThreads,threads,data2);
        my_join(threads,numThreads);
        //printf("pa\n");
        //PointMatrix(pa);

        pap=vectorMul(pa,p1);
        //printf("pap=%f\n",pap);

        al=r1r/pap;
        //printf("alpha=%f\n",al);

        mulReal(al,p1,ap1);
        //printf("ap1\n");
        //PointMatrix(ap1);

        addMatrix(X,ap1,x2);
        //printf("x2\n");
        //PointMatrix(x2);

        mulReal(al,A,aA);
        //printf("aA\n");
        //PointMatrix(aA);

        initData(numThreads,data3,aA,p1,aAp);
        my_thread(numThreads,threads,data3);
        my_join(threads,numThreads);
        //printf("aAp\n");
        //PointMatrix(aAp);

        subMatrix(r1,aAp,r2);
        //printf("r2\n");
        //PointMatrix(r2);

        fr=finish(r2);
        //printf("fr=%f\n",fr);

        if(fr<0.001)
            break;
        
        transpose(r2,r2t);
        //printf("r2t\n");
        //PointMatrix(r2t);
        
        r2r=vectorMul(r2t,r2);
        //printf("r2r=%f\n",r2r);

        bat=r2r/r1r;
        //printf("bat=%f\n",bat);
        mulReal(bat,p1,bp);
        addMatrix(r2,bp,p2);


        fuzhi(p2,p1);
        fuzhi(r2,r1);
        fuzhi(x2,X);

        //printf("p1\n");
        //PointMatrix(p1);
        //printf("r1\n");
        //PointMatrix(r1);
        //printf("X\n");
        //PointMatrix(X);

        //printf("迭代分界线\n");
        i++;
        //PointMatrix(x2);
    }

    gettimeofday(&t2,NULL);
    printf("Количество итерации:%d\n",i);
    double Xn=B.m[B.x*B.y-1]/A.m[A.x*A.y-1];
    double Xn_1=(B.m[B.x*B.y-2]-A.m[(A.x-1)*A.y-1]*Xn)/(A.m[(A.x-1)*A.y-2]);
    printf("Метод Гаусса:x(n):%f,x(n-1)=%f\n",Xn,Xn_1);
    printf("Метод сопряженных градиентов:x(n)=%f,x(n-1)=%f\n",x2.m[x2.x*x2.y-1],x2.m[x2.x*x2.y-2]);
    double time=t2.tv_sec-t1.tv_sec+(double)(t2.tv_usec-t1.tv_usec)/1000000;
    printf("Размерность матрицы:%d\n",dim);
    printf("Количество потоков:%d\nВремя:%fs\n",numThreads,time);
    

    pthread_exit(NULL);

}