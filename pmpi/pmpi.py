import numpy as np
import datetime
import math
from mpi4py import MPI

dim = 1000

class Matrix:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.m = np.zeros((x, y))

    def positive_definite_matrix(self):
        sum = 1.0
        for i in range(0, self.x):
            for j in range(0, self.y):
                if j < i:
                    self.m[i][j] = 0
                else:
                    self.m[i][j] = sum
                    sum += 1

    def transpose(self, t):
        for i in range(0, self.x):
            for j in range(0, self.y):
                self.m[i][j] = t.m[j][i]

    def vector_mul(self, b):
        sum = 0.0
        for i in range(0, self.y):
            sum += self.m[0][i] * b.m[i][0]
        return sum

    def mul_real(self, r, rm):
        for i in range(0, self.x):
            for j in range(0, self.y):
                self.m[i][j] = r * rm.m[i][j]

    def add_matrix(self, m1, m2):
        for i in range(0, self.x):
            for j in range(0, self.y):
                self.m[i][j] = m1.m[i][j] + m2.m[i][j]

    def sub_matrix(self, m1, m2):
        for i in range(0, self.x):
            for j in range(0, self.y):
                self.m[i][j] = m1.m[i][j] - m2.m[i][j]

    def fuzhi(self, M):
        for i in range(0, self.x):
            for j in range(0, self.y):
                self.m[i][j] = M.m[i][j]

    def finish(self):
        lin = 0.0
        sum = 0.0
        for i in range(0, self.x):
            for j in range(0, self.y):
                lin = abs(self.m[i][j])
                sum += lin * lin
        return sum

    def resetmatrix(self):
        for i in range(0, self.x):
            for j in range(0, self.y):
                self.m[i][j] = 0

    def mult_matrix(self, m1, m2):
        for i in range(0, self.x):
            for j in range(0, self.y):
                sum = 0.0
                for k in range(0, m1.y):
                    sum += m1.m[i][k] * m2.m[k][j]
                self.m[i][j] = sum

    def reset(self):
        for i in range(0, self.x):
            for j in range(0, self.y):
                self.m[i][j] = 0

def MultMatrix(id,numprocs,mul1,mul2,prod):
    prod.reset()
    num=math.ceil(dim/numprocs);
    indexProd=id*num;
    for i in range(0,num):
        if(indexProd+i==prod.x*prod.y):
            break
        X=int((indexProd+i)/prod.y);
        Y=(indexProd+i)%prod.y;
        sum=0.;
        for j in range(0,mul1.y):
            sum+=mul1.m[X][j]*mul2.m[j][Y];
        prod.m[X][Y]=sum;
    if (id != 0):
        comm.send(prod.m,dest=0, tag=99)
    else:
        for k in range(1,numprocs):
            m=comm.recv(source=k, tag=99)
            for i in range(0,prod.x):
                for j in range(0,prod.y):
                    prod.m[i][j]+=m[i][j]



t1=datetime.datetime.now()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

A=Matrix(dim,dim);
A.positive_definite_matrix()
X=Matrix(dim,1);
B=Matrix(dim,1);
for i in range(0,dim):
    B.m[i][0]=2

Ax = Matrix(dim, 1);
r1 = Matrix(dim, 1);
rt = Matrix(1, dim);
p1 = Matrix(dim, 1);
pt = Matrix(1, dim);
pa = Matrix(1, dim);
ap1 = Matrix(dim, 1);
x2 = Matrix(dim, 1);
aA = Matrix(dim, dim);
aAp = Matrix(dim, 1);
r2 = Matrix(dim, 1);
r2t = Matrix(1, dim);
bp = Matrix(dim, 1);
p2 = Matrix(dim, 1);

Ax.mult_matrix(A,X)
r1.sub_matrix(B,Ax)
p1.fuzhi(r1)
i=0
while i<20:
    rt.transpose(r1)
    r1r=rt.vector_mul(r1)
    pt.transpose(p1)

    MultMatrix(rank,size,pt,A,pa)
    pa.m=comm.bcast(pa.m, root=0)
    comm.barrier() 

    pap=pa.vector_mul(p1)
    al=r1r/pap

    ap1.mul_real(al,p1)

    x2.add_matrix(X,ap1)
    aA.mul_real(al,A)

    MultMatrix(rank,size,aA,p1,aAp)
    aAp.m=comm.bcast(aAp.m, root=0)
    comm.barrier() 
    
    r2.sub_matrix(r1,aAp)
    fr=r2.finish()
    if fr<0.01:
        break

    r2t.transpose(r2)

    r2r=r2.m[0][0]*r2.m[0][0]

    bat=r2r/r1r
    bp.mul_real(bat,p1)
    p2.add_matrix(r2,bp)
    p1.fuzhi(p2)
    r1.fuzhi(r2)
    X.fuzhi(x2)
    i+=1

comm.barrier()
t2=datetime.datetime.now()
#print(X.m)
if(rank==0):
    #Xn=B.m[B.x-1][B.y-1]/A.m[A.x-1][A.y-1];
    #Xn_1=(B.m[B.x-1][B.y-2]-Xn*A.m[A.x-2][A.y-1])/A.m[A.x-2][A.y-2]
    #print("Метод Гаусса:x(n):%f,x(n-1)=%f\n"%(Xn,Xn_1));
    #print("Метод сопряженных градиентов:x(n)=%f,x(n-1)=%f\n"%(x2.m[x2.x-1][x2.y-1],x2.m[x2.x-1][x2.y-2]));
    print("Количество итерации%d"%i)
    time=(t2-t1).microseconds/1000000
    print("Размерность матрицы:%d\n"%dim);
    print("Количество потоков:%d\nВремя:%fs\n"%(size,time));







