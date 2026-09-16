#Ax = b 
N = 20
delta = 0.05

A = operator_A(N)
b = np.zeros(N)
for i in range(N):
    sumA = 0
    for j in range(N):
        sumA += A[i][j]
    b[i] = sumA 

x = np.linspace(0,1,N)
sigma = np.random.normal(N)

b_delta = b + delta * sigma

y = np.linalg.solve(A, b_delta)
plt.plot(y)
