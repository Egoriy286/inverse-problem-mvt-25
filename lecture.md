# 1

ОДУ паробалического типа

$(1) \quad \dfrac{\partial u}{\partial t} = \dfrac{\partial^2 u}{\partial x^2} + p(t)u, \quad x \in (0,l), \, t \in (0, \overline{t}]$

$(2) \quad u(0,t) = \mu_1(t), \quad u(l,t) = \mu_2(t), \quad t \in (0, \overline{t}]$

$(3) \quad u(x,0) = u_0(x), \quad x\in[0,l]$

Ставим обратную задачу 

ставим условие переопределения

$(4) \quad u(\overline{x}, t) = \phi(t), \quad \overline{x} \in (0,l), t \in (0, t] , q(0) = u_0(\overline{x})$


![alt text](<graph 1.png>)


**Ищем** 

$ u(x,t)-? $

чтобы найти $p(t)$ мы задаем условие переопределения

# 2


Дискретизация

$ x_i = ih, \quad i = \overline{0, N}, \quad h = \frac{l}{N}$

$ \overline{x} = nh$

$ t_j = j \tau, \quad j = \overline{0, J}, \tau = \frac{\overline{t}}{J} $


Конечно-разностная апроксимация

Кранка Николсона 

$(5) \quad \dfrac{y_{i}^{j} - y_{i}^{j+1}}{\tau} = \dfrac{y_{i+1}^{j+1} - 2y_{i}^{j+1} + y_{i-1}^{j+1}}{2h^2} + \dfrac{y_{i+1}^{j} - 2y_{i}^{j} + y_{i-1}^{j}}{2h^2} - \dfrac{1}{2} \left(p^{j+1}y_{i}^{j} + p^{j}y_{i}^{j+1} \right) \\ y_{0}^{j} = \mu_1^{j}, \quad y_{N}^{j}= \mu_{2}^{j} $

$(6) \quad y_{i}^{0} = u_0 (x_i), \quad i = \overline{0, N}$

$(7) \quad y_n^j = \phi^j, \quad j = \overline{1, J}, \phi^0 = u_0(x_n)$



$i=\overline{1, N-1}, \quad j = \overline{1, J};$

Точность апроксимации:

$\psi = O(\tau^2 + h^2)$

Решение (5) ищем в виде линейной формы:

$(8) \quad y_i^j = z_i + p^j v_i, \quad i = \overline{0,N}$

coefs for tridiagonal

$ r = \dfrac{2h^2}{\tau}$

$a_i = 1 = b_i$

$c_i = 2 + r + p^{j} h^2 p^{j-1}$

$d_i = ry_{i}^{j-1} + y_{j+1}^{j-1} - 2y_{j}^{j-1} + y_{j-1}^{j-1} + h^2 p^{j-1} y_i^{j}$

$\overline{d_i} = ry_i^{j-1} + y_{i+1}^{j-1}-2y_{i}^{j-1} + y_{j-1}^{j-1}$

$\overline{c_i} = 2 + r$ 

$(9) \quad y_{0}^{j} = \mu_1^{j}, \quad  a_i y_{i-1}^{j} - c_{i}y_{i}^{j} + b_i y_{i-1}^{j} + d_i = 0 , \quad i = \overline{1, N-1}, y_N^{j} = \mu_2^{j} $

Итого:

$ \{ a_i z_{i-1} - \overline{c_i} z_i + b_i z_{i-1} + \overline{d_i} \} + p^{j} \{ a_i v_{i-1} - (\overline{c_i} + h^2 p^{j-1})v_{i} + b_{i} v_{i-1} + h^2y_i^{j-1} \} = 0 , \quad i = \overline{1, N-1}$


$(10) \quad z_0 = \mu_1, \quad a_i z_{i-1} - \overline{c_i} z_{i} + b_{i} z_{i+1} + d_i = 0, \quad i = \overline{1, N-1} ,\quad z_N = \mu_2 $

$(11) \quad v_0 = 0, \quad a_i v_{i-1} + \overline{c_i} v_i + b_i v_{i+1} + h^2 y_{i}^{j-1} = 0, \quad i = \overline{1, N-1} , v_N = 0 $

$(12) \quad y_{N}^{i} \to z_N + p^j v_n = \phi_j \to p^j = \dfrac{\phi_j - z_N}{v_N}$

$(6) \to (10) \to (11) \to (12)$