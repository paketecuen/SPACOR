
## Case 2: Parallel RL load switching the resitive part with Sinusoidal source

Given a parallel RL circuit with a switch placed before of the **R**, it is supplied with any voltage giving rise to a current. Through some electric meter, you have access to both v and i. It is desired to estimate the value of R and L for every time instant, assuming that these values can change over time. The interval of switching are <img src="https://render.githubusercontent.com/render/math?math=0 < t < 2"> switch off, <img src="https://render.githubusercontent.com/render/math?math=2 < t < 6">  switch on,  <img src="https://render.githubusercontent.com/render/math?math=6 < t < 13"> switch off,  <img src="https://render.githubusercontent.com/render/math?math=13 < t < 18"> switch and, finally, <img src="https://render.githubusercontent.com/render/math?math=18 < t < 20"> switch off.



Applying the Kirchhoff's first law, or KCL,
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=i = i_G + i_L = G_v + \frac{1}{L}\int v = Gv +\Gamma \check{v}"></p>


where the symbol has been used to indicate the voltage integral and, in addition, the symbol  for the inverse of the inductance. In the above expression, it has been assumed (for simplicity) that the initial conditions are such that they give rise to null values.

The use of the normal and tangent vector define the concept of the osculating plane. The variation of this plane in space is what allows us to correctly estimate the parameters R and L. . For the 2D case at hand, it is sufficient to use the vector y (since it is contained in the osculating plane itself, although in 3D it will be necessary to use the normal vector since we are not sure that y is included in the osculating plane) and the tangent vector, calculated as the derivative of y. Thus, we have,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=y = v\sigma_1 + \check{v}\sigma_2 + i \sigma_3"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=y' = v'\sigma_1 + v\sigma_2 + i' \sigma_3"></p>

and proceed to calculate the oscillating plane,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=K_{osc} = y \wedge y' = (v^2 - \check{v}v')\sigma_{12} + (\check{v}i'-iv)\sigma_{23} + (iv'-vi')\sigma_{31}"></p>

Following the reasoning of finding the ratio between <img src="https://render.githubusercontent.com/render/math?math=\sigma{13}">  and <img src="https://render.githubusercontent.com/render/math?math=\sigma{12}"> , as well as <img src="https://render.githubusercontent.com/render/math?math=\sigma{31}">  and <img src="https://render.githubusercontent.com/render/math?math=\sigma{12}">  the value of <img src="https://render.githubusercontent.com/render/math?math=\G">  and <img src="https://render.githubusercontent.com/render/math?math=\Gamma">  is now,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=G = \frac{iv-\check{v}i'}{v^2 -\check{v}v'}"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=\Gamma = \frac{iv-\check{v}i'}{v^2 -\check{v}v'}"></p>

[Matlab simulation for case 2](examples/matlabCase2.md)
