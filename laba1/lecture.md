# 1. Постановка обратной задачи

Уравнение параболического типа:

$(1)\quad \dfrac{\partial u}{\partial t} = \dfrac{\partial^2 u}{\partial x^2} + p(t)\,u,\quad x\in(0,l),\ t\in(0,T]$

Граничные условия:

$(2)\quad u(0,t)=\mu_1(t),\quad u(l,t)=\mu_2(t),\quad t\in(0,T]$

Начальное условие:

$(3)\quad u(x,0)=u_0(x),\quad x\in[0,l]$

Условие переопределения (внутреннее измерение):

$(4)\quad u(\bar{x},t)=\phi(t),\quad \bar{x}\in(0,l),\ t\in[0,T],\ \phi(0)=u_0(\bar{x})$

Требуется найти одновременно:

$u(x,t)$ и $p(t)$.

![alt text](<graph 1.png>)

# 2. Дискретизация и разностная схема

Сетка:

$x_i=ih,\ i=\overline{0,N},\ h=\dfrac{l}{N},\ \bar{x}=x_n=nh$

$t_j=j\tau,\ j=\overline{0,J},\ \tau=\dfrac{T}{J}$

Обозначения:

$y_i^j \approx u(x_i,t_j),\quad p^j \approx p(t_j),\quad \mu_k^j=\mu_k(t_j),\quad \phi^j=\phi(t_j)$

Схема Кранка-Николсона с линейризацией члена $p(t)u$:

$$(5)\quad
\dfrac{y_i^j-y_i^{j-1}}{\tau}
=
\dfrac{1}{2h^2}
\Big[(y_{i+1}^j-2y_i^j+y_{i-1}^j)+(y_{i+1}^{j-1}-2y_i^{j-1}+y_{i-1}^{j-1})\Big]
+\dfrac{1}{2}\Big(p^j y_i^{j-1}+p^{j-1}y_i^j\Big),$$

$i=\overline{1,N-1},\ j=\overline{1,J}.$

Граничные, начальные и дополнительные условия в сеточной форме:

$(6)\quad y_0^j=\mu_1^j,\quad y_N^j=\mu_2^j,\quad j=\overline{1,J}$

$(7)\quad y_i^0=u_0(x_i),\quad i=\overline{0,N}$

$(8)\quad y_n^j=\phi^j,\quad j=\overline{1,J},\quad \phi^0=u_0(x_n)$

Порядок аппроксимации: $O(\tau^2+h^2)$.

# 3. Трёхдиагональная форма

Положим:

$r=\dfrac{2h^2}{\tau},\quad a_i=b_i=1,\quad c_i=2+r-h^2p^{j-1}$

$\bar d_i=r\,y_i^{j-1}+(y_{i+1}^{j-1}-2y_i^{j-1}+y_{i-1}^{j-1})$

Тогда (5) эквивалентно:

$(9)\quad a_i y_{i-1}^j-c_i y_i^j+b_i y_{i+1}^j+\bar d_i+h^2p^j y_i^{j-1}=0,\quad i=\overline{1,N-1}.$

# 4. Разложение по неизвестному $p^j$

Ищем $y_i^j$ в виде:

$(10)\quad y_i^j=z_i^j+p^j v_i^j,\quad i=\overline{0,N}.$

После подстановки (10) в (9) получаем две линейные задачи:

Для $z^j$:

$(11)\quad
\begin{cases}
z_0^j=\mu_1^j,\\
a_i z_{i-1}^j-c_i z_i^j+b_i z_{i+1}^j+\bar d_i=0,\quad i=\overline{1,N-1},\\
z_N^j=\mu_2^j.
\end{cases}$

Для $v^j$:

$(12)\quad
\begin{cases}
v_0^j=0,\\
a_i v_{i-1}^j-c_i v_i^j+b_i v_{i+1}^j+h^2y_i^{j-1}=0,\quad i=\overline{1,N-1},\\
v_N^j=0.
\end{cases}$

Из условия переопределения в узле $n$:

$(13)\quad y_n^j=z_n^j+p^j v_n^j=\phi^j
\ \Longrightarrow\ 
p^j=\dfrac{\phi^j-z_n^j}{v_n^j},\quad v_n^j\neq 0.$

Далее:

$(14)\quad y_i^j=z_i^j+p^j v_i^j,\quad i=\overline{0,N}.$

# 5. Алгоритм по времени

$(7)\rightarrow(11)\rightarrow(12)\rightarrow(13)\rightarrow(14)$ для каждого $j=1,\dots,J$.
