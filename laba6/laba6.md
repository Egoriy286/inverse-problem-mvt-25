Решить с помощью метода конечных элементов $u( \vec{x},t) \quad p(t) ?$

$$(1) \quad  \dfrac{\partial u}{\partial t} + Lu =f \quad \vec{x} \in \Omega , \quad t \in (0, T]$$

$$(2) \quad Lu = - \Delta u - p u  $$

$$(3) \quad u|_{\Gamma} = 0, \quad \vec{x} \in \partial \Omega, \quad t \in (0, T] $$

$$(4) \quad u(\vec{x}, 0) = u_0(\vec{x}), \quad \vec{x} \in \Omega $$

$$(5) \quad p(t): \int_{\Omega} \mu(t) \cdot u(\vec{x}, t) = \phi(t) , \quad t \in [0,T]$$

---
обратная разностная схема 

$$(6) \quad \dfrac{y-\check{y}}{\tau} + Ly - p(t)\cdot \check{y} = f , \quad t = t_j , \quad \vec{x} \in \Omega$$

$$(7) \quad y = z +p(t)\cdot \omega , \quad \vec{x} \in {\Omega} , t = t_j $$

подставляем (7) на уравнение (6) получаем

$$\dfrac{z - \check{y}}{\tau} + p(t) \dfrac{\omega}{\tau} - L z - pL\omega - p \check{y} = f$$

$\check{y} - у на нижем слое$

$$ \left( \dfrac{z - \check{y}}{\tau} - Lz = f \right) + p \left( \dfrac{\omega}{\tau} - L\omega - \check{y} \right) = 0 \Rightarrow  $$

$$ \Rightarrow -z + \tau L z + f + \check{y} = 0 \qquad(8) \\
    - \omega + \tau L \omega - y = 0,  \qquad\qquad (9)
$$

---
шаг действий

$$ (10) \Rightarrow (8) - (9) $$
$$y^{0} \to z,w, \quad \vec{x} \in \overline{\omega}$$