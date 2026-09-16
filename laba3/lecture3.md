$$\begin{align} p \dfrac{\partial u}{\partial t} = \dfrac{\partial^2 u}{\partial x^2},\quad x \in (0, l) ,\quad t \in (0, T]  \end{align}$$

$$u(x,t), p(t) - ? \quad r = \dfrac{h^2}{\tau}, \quad c= 2+r$$

$$\begin{align} u(0,t) = \mu_1(t), \quad u(l,t) = \mu_2(t),\quad t \in (0, T] \end{align}$$

$$\begin{align} u(x,0) = u_0(x) , \quad x \in [0, l] \end{align}$$

$$\begin{align} \int_{0}^{l} u(x,t)dx = \phi , \quad t \in [0,T] \end{align}$$

Решение обратной задачи на квази решение
(3) => (5)

$$\begin{align*} y=u_0(x_i), \quad i = \overline{0,n} , \phi^{0} = \dfrac{h}{2} \sum_{i=0}^{n} y_{i}^{0} \end{align*}$$

$$\begin{align*} u(0,t) = \mu_1(t), \quad u(l,t) = \mu_2(t),\quad t \in (0, T]  \end{align*}$$


----
квазирешение

$$(3)' \to  y_{i}^{0} = u_0(x_i), \quad i = \overline{0,n}$$


$$(1-2)' \to  y_{j}^{0} = \mu_{1}^{j}, \quad p^{j} \dfrac{y^j - \hat{y}}{\tau} =  \dfrac{y_{i+1} - 2 y_{i} + y_{i-1}}{h^2}, \quad i = \overline{1,n-1}; \quad y_n = \mu_2; \quad j = \overline{1, M}$$


находим квазирешение
$$j = \overline{0, J-1}$$
$$\begin{align} y_{i}^{0} = u_{0}(x_i), \quad i = \overline{0,n}, \quad \phi^0 = \dfrac{h}{2} \sum_{i=1}^{n} y_i^0 \end{align}$$


$$\begin{align} \begin{cases} y_0 = \mu_1, \quad y_{i-1} - (2+ Pr)y_i + y_{i+1}+ P r \hat{y_i} = 0, \quad i=\overline{1, n-1} ; \quad  y_n = \mu_2, \\
\phi^{j} = \dfrac{h}{2} \sum_{i=1}^{n}(y_{i-1} + y_{i}) \end{cases} \end{align}$$



интеграл (4) найдем дискретный аналог при заданном p()

формула трапеций, правых и левых прямоугольников будет совпадать


мы не знаем что такое phi_i из прямой задачи опрделяем условие опредление 


---
переходим к обратной задаче


$$ j=\overline{1, J} $$



$$\begin{align} \begin{cases} y_{0} = \mu_1, \quad y_{i-1} - (2 + rp)y_i + y_{i+1} + rp\hat{y}, \quad i = \overline{1, n-1}; \quad y_n= \mu_2, \\
 \dfrac{h}{2} \sum_{i=1}^{n}(y_{i-1} + y_{i}) = \phi \end{cases} \end{align}$$


(10)

$$\begin{align} y_i = z_i + p \cdot w_i , \quad i = \overline{0,n} \end{align}$$

$$ y_0 = \mu_1 \to w_0 , z_0 = \mu_1 ; \quad y_n = \mu_2 \to w_n =0, z_n = \mu_2$$


примерно решение погрешность $10^{-12}$


---



$$ p \dfrac{\partial u}{\partial t} = \dfrac{\partial p u}{ \partial t} - u \dfrac{\partial p}{\partial t} $$


$$ p \dfrac{\partial u}{\partial t} + u \dfrac{\partial p}{\partial t} = \dfrac{\partial p u}{ \partial t}  $$

$$ \dfrac{pu - \hat{p} u}{\tau} = \hat{p}\dfrac{u - \hat{u}}{\tau} + \hat{u} \dfrac{p - \hat{p}}{\tau} = \dfrac{\hat{p}u - \hat{p} \hat{u} + \hat{u} p - \hat{u} \hat{p}}{\tau} \to$$

$$ \hat{p} \dfrac{u - \hat{u}}{\tau} = \dfrac{p u - \hat{p} u}{\tau} - \dfrac{p \hat{u} - p \hat{u}}{\tau} $$

$$ ((z_{i-1} - (2 - r\hat{p}) z_i + z_{i+1} + r \hat{p} \hat{y}) + p(w_{i-1} + (2+r\hat{p}) w_i + w_{i-1})  $$