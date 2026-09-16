import numpy as np; import matplotlib.pyplot as plt                                     # Студента Монастырева Егория ПМИ-21 4-курс  
v = 1.0; D = 1.e-3; l = 1; Time = 0.8; N=500; h = l/N                                   # Исходные данные ипараметры сетки 
#if(tau<= st): print("Устойчива")
#else: print("Неустойчива")
M=400; tau=Time/M; r=D*tau/(h*h); e=v*tau/h; st = h**2/(2.*D*tau); print(tau,h,N,M,r,st)
x=np.zeros(N+1,'float'); T=np.zeros((M+1,N+1),'float'); y = np.zeros(N+1,'float')
for i in range(0,N+1): T[0,i] = 0.; x[i]=i*h                                            # Начальное условие. Расщепление по физическим процессам
for j in range(1, M+1):                                                                       
    T[j,0] = 1.; T[j,N] = 0.; y[0]=1.                                                   # Граничные условия
    for i in range(1,N+1): y[i] = T[j-1,i-1]
    for i in range(1, N): T[j,i] = .5*(y[i-1]+y[i+1])

M=800; tau=Time/M; r=D*tau/(h*h); e=v*tau/h; st = h**2/(2.*D*tau); print(tau,h,N,M,r,st); T1= np.zeros((M+1, N+1))  # Временных слоее 2 раза больше
for i in range(0,N+1): T1[0,i] = 0.                                                     # Явная схема. Начальное условие
for j in range(1, M+1): 
    T1[j,0] = 1.; T1[j,N] = 0.; y[0] = 1.                                               # Граничные условия
    for i in range(1, N+1): y[i]=T1[j-1,i] 
    for i in range(1,N): T1[j,i]=y[i]+r*(y[i+1]-2*y[i]+y[i-1])-e*(y[i]-y[i-1])
    
fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12, 4))
im1 = ax1.imshow(T, extent=[0, l, 0, Time], aspect='auto', cmap='hot', origin='lower', interpolation='nearest')
ax1.set_ylabel("Время (t)"); ax1.set_xlabel("Координата (x)"); ax1.set_title("По физическим процессам"); ax1.grid(True)
fig.colorbar(im1, ax=ax1, orientation='vertical', label='U(x,t)')

im2 = ax2.imshow(T1, extent=[ 0, l,0, Time], aspect='auto', cmap='hot', origin='lower', interpolation='nearest')
ax2.set_ylabel("Время (t)"); ax2.set_xlabel("Координата (x)"); ax2.set_title("Явная схема")
fig.colorbar(im2, ax=ax2, orientation='vertical', label='U(x,t)'); ax2.grid(True); plt.show()

'''fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4)); 
for j in range(0, M, M // 10): ax1.plot(x, T[j, :], label=f't={j * tau:.2f}')
ax1.set_ylabel("U(x,t)"); ax1.set_xlabel("Координата (x)"); ax1.set_title("По физическим процессам"); ax1.legend(); ax1.grid(True)

for j in range(0, M, M // 10): ax2.plot(x, T1[j, :], label=f't={j * tau:.2f}')
ax2.set_ylabel("U(x,t)"); ax2.set_xlabel("Координата (x)"); ax2.set_title("Явная схема"); ax2.grid(True); ax2.legend(); plt.show()'''