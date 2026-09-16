## Обратная задача для уравнения переноса-диффузии

Найти функции $u(x,t)$ и $v(t)$, такие что

$$
\frac{\partial u}{\partial t}
=
\frac{\partial^2 u}{\partial x^2}
- v(t)\frac{\partial u}{\partial x},
\qquad
x\in(0,l),\; t\in(0,\overline{t}],
$$

$$
u(0,t)=\mu_1(t),\qquad u(l,t)=\mu_2(t),
\qquad t\in(0,\overline{t}],
$$

$$
u(x,0)=u_0(x),\qquad x\in[0,l],
$$

$$
\int_0^l u(x,t)\,dx=\phi(t),
\qquad t\in(0,\overline{t}],
$$

$$
v(t)>0.
$$

### Сетка и обозначения

$$
x_i=ih,\; i=\overline{0,N},\qquad
t_j=j\tau,\; j=\overline{0,J},
$$

$$
h=\frac{l}{N},\qquad
\tau=\frac{\overline{t}}{J},\qquad
y_i^j\approx u(x_i,t_j),\qquad
y_i^0=u_0(x_i).
$$

Вводим константы:

$$
r=\frac{2h^2}{\tau},\qquad
c=2+r,\qquad
e=2-r.
$$

### Разностная схема (прямая задача)

Для $j=\overline{0,J-1}$, $i=\overline{1,N-1}$:

$$
y_{i+1}^{j+1}-c\,y_i^{j+1}+y_{i-1}^{j+1}
=
-y_{i+1}^{j}+e\,y_i^{j}-y_{i-1}^{j}
+2h\,v^{j+1}\bigl(y_i^{j}-y_{i-1}^{j}\bigr),
$$

$$
y_0^{j+1}=\mu_1^{j+1},\qquad
y_N^{j+1}=\mu_2^{j+1}.
$$

### Линеаризация для обратной задачи

Представление:

$$
y_i^{j+1}=z_i^{j+1}+v^{j+1}w_i^{j+1}.
$$

Тогда получаем две трёхдиагональные задачи.

Для $z$:

$$
z_{i+1}^{j+1}-c\,z_i^{j+1}+z_{i-1}^{j+1}
=
-y_{i+1}^{j}+e\,y_i^{j}-y_{i-1}^{j},
\qquad i=\overline{1,N-1},
$$

$$
z_0^{j+1}=\mu_1^{j+1},\qquad
z_N^{j+1}=\mu_2^{j+1}.
$$

Для $w$:

$$
w_{i+1}^{j+1}-c\,w_i^{j+1}+w_{i-1}^{j+1}
-2h\bigl(y_i^{j}-y_{i-1}^{j}\bigr)=0,
\qquad i=\overline{1,N-1},
$$

$$
w_0^{j+1}=0,\qquad
w_N^{j+1}=0.
$$

### Восстановление $v^{j+1}$ из условия переопределения

$$
\phi^{j+1}
=
\sum_{i=1}^{N}\frac{h}{2}\left(y_i^{j+1}+y_{i-1}^{j+1}\right).
$$

Подставляя $y=z+v^{j+1}w$, получаем

$$
v^{j+1}
=
\frac{
\phi^{j+1}
-\sum_{i=1}^{N}\frac{h}{2}\left(z_i^{j+1}+z_{i-1}^{j+1}\right)
}{
\sum_{i=1}^{N}\frac{h}{2}\left(w_i^{j+1}+w_{i-1}^{j+1}\right)
}.
$$

Если знаменатель близок к нулю, на практике берут $v^{j+1}=v^j$.

### Оценка начального значения $v^0$

$$
\phi_t(0)=u_x(l,0)-u_x(0,0)-v^0\bigl(u(l,0)-u(0,0)\bigr),
$$

$$
v^0=
\frac{u_x(l,0)-u_x(0,0)-\phi_t(0)}
{u(l,0)-u(0,0)}.
$$

### Последовательность шага по времени

$$
(z,w)\;\rightarrow\;v^{j+1}\;\rightarrow\;y^{j+1}.
$$

Источник: Безношенко, *Сибирский журнал вычислительной математики и математического моделирования*.
