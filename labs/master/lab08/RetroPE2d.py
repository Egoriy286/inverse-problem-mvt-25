import numpy as np; from LibrRetro2d import u0,poisson,gr,gr1                   #Iterative identification of u0(x). D/S winh sigma 2d
from numpy import linalg as LA; from numpy import matrix
from numpy import linalg; import matplotlib.pyplot as plt
N=50; M=40; l=16.; h=l/N; t0=.25; tt=1.; tau=(tt-t0)/M; si=.5-h*h/(12.*tau)
r1=tau*si/(h*h); r2=tau*(1.-si)/(h*h); kk=190; eps=1.e-12
rf=np.zeros((kk),'float');       tr=np.zeros((kk),'float'); rf[0]=0.
f=np.zeros((N+1)*(N+1),'float'); y=np.zeros((N+1)*(N+1),'float')
z=np.zeros((N+1)*(N+1),'float'); p=np.zeros((N+1)*(N+1),'float')
v=np.zeros((N+1)*(N+1),'float'); r=np.zeros((N+1)*(N+1),'float')
u=np.zeros((N+1)*(N+1),'float'); w=np.zeros((N+1)*(N+1),'float')

A=poisson(N+1,-r1,1.+4.*r1); B=A.toarray(); C1=linalg.inv(B)
A1=poisson(N+1,r2,1.-4.*r2); B1=A1.toarray(); C=np.dot(C1,B1)
for j in range(0,N+1):                                                          #Exact solution
    for i in range(0,N+1): k=i+(N+1)*j; f[k]=u0((i*h-.5*l),(j*h-.5*l),t0); w[k]=u0((i*h-.5*l),(j*h-.5*l),tt); y[k]=f[k]
for m in range(1,M+1): y=np.dot(C,y)
e=np.dot(f,f); Y=np.reshape(y,(N+1,N+1)); W=np.reshape(w,(N+1,N+1)); gr1(N,l,W,1,1)
F=np.reshape(f,(N+1,N+1)); gr1(N,l,F,1,1); gr1(N,l,W-Y,1,1); v=w; u=w

for m in range(1,M+1): u=np.dot(C,u)                                            #Iterative method CG
r=w-u; p=r; a=np.dot(r,r)

for K in range(1,kk):                                                           #Iterations
    z=p
    for m in range(1,M+1): z=np.dot(C,z)
    c=np.dot(z,p); a1=a; alp=a1/c; v=v+alp*p; r=r-alp*z; a=np.dot(r,r)
    if np.sqrt(a)<eps: break
    p=r+a*p/a1; cc=np.sqrt(np.dot(f-v,f-v)/e); rf[K]=cc; print (K, cc,alp, a/a1)

print (K, N, M, tt)
V=np.reshape(v,(N+1,N+1)); gr1(N,l,V,1,1); gr1(N,l,V-F,1,1); gr(N,l,-.8e-2,.4e-2,21,F-V)