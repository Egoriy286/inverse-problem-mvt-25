? $u(x,t), \quad v(t):\infty > v(t) > 0$
<!-- $$\begin{align}  \end{align}$$ -->
$$\begin{align} \dfrac{\partial u}{\partial t} = \dfrac{\partial^2 u}{\partial x^2} - v \dfrac{\partial u}{\partial x}, \qquad x \in (0, l), \quad t \in (0,\overline{t}] \end{align}$$

$$\begin{align} u(0, t) = \mu_1(t) , \quad u(l,t) = \mu_2(t), \quad t \in (0, \overline{t}] \end{align}$$

$$\begin{align} u(x,0) = u_0(x) , \quad x \in [0, l] \end{align}$$

$$\begin{align} \int_{0}^{l} u(x,t)dx = \phi(t), t \in (0, \overline{t}]  \end{align}$$

$$\begin{align} y_{i}^{0} = u_0(x_i), \quad i = \overline{0, N} \end{align}$$


определение констант:
$$r = \dfrac{2h^2}{\tau}, \quad c = 2 + r, \quad e = 2 -r, \quad g = vh$$

МКР апроксимация Кранка-Николсон (черновик):

$$r (y_{i}^{j} - y_{i}^{j-1}) = \Big[(y_{i+1}^j-2y_i^j+y_{i-1}^j)+(y_{i+1}^{j-1}-2y_i^{j-1}+y_{i-1}^{j-1})\Big] \\ - g (y_{i+1}^{j} + y_{i+1} - y_{i-1}^{j} - y_{i-1}^{j} - y_{i-1}^{j}) \cdot(v_{i}^{j} + v_{i}^{j-1})$$

будет квадрат при подставки (7) меняем апроксимацию

$$r (y_{i}^{j} - y_{i}^{j-1}) = \Big[(y_{i+1}^j-2y_i^j+y_{i-1}^j)+(y_{i+1}^{j-1}-2y_i^{j-1}+y_{i-1}^{j-1})\Big] \\ - v^j 2h\left(y^{j-1}-y_{i-1}^{j-1}\right)$$




$y_0^j = \mu_1^{j}$

$y_N^{j} = \mu_2^{j}$

$$\begin{align} y_{i+1}^{j} - c y_{i}^{j} + y_{i-1}^{j} + y_{i+1}^{j-1} - e y_{i}^{j-1}  y_{i-1}^{j-1} - 2v^{j} h(y_i^{j-1} - y_{i-1}^{j-1})\end{align}$$
проходимся $i=\overline{1, N-1}$


переходим к обратной апроксимации

$$\begin{align} y_{i}^{j} = z_i + v^j \cdot w_i\end{align}$$

$$\begin{align} z_0^{j} = \mu_1^j, \quad z_{i+1} - cz_i + z_{i-1} = - y_{i+1}^{j-1} + e y_i^{j-1} - y_{i-1}^{j-1} \end{align}$$
$i = \overline{1, N-1}, \quad z_N = \mu_2^i$

$$\begin{align} w_0^j = 0, \quad w_{i+1} - c w_i + w_{i-1} - 2h (y_i^{j-1} - y_{i-1}^{j-1}) = 0 \end{align}$$
$i = \overline{1, N-1}, \quad w_N = \mu_2^j$

$$\begin{align} \sum_{i=1}^{N} \dfrac{h}{2} (y_i + y_{i-1}) = \phi^{j} \to v_j = \dfrac{\phi^j - \sum_{i=1}^{N} \dfrac{h}{2} (z_i + z_{i-1}) }{\sum_{i=1}^{N} \dfrac{h}{2} (w_i + w_{i-1}) } \end{align}$$

правило
$(5) \to (8,9) \to (10) \to (7)$

Безношенко Сибирский журнал матмоделирования