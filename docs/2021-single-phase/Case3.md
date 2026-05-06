## Case 3: Serie RL with Sinusoidal source

Given a parallel RL circuit as in Figure 3, it is supplied with any voltage giving rise to a current. Through some electric meter, you have access to both v and i. It is desired to estimate the value of **R** and **L**for every time instant, assuming that these values can change over time. Similarly, the method can also be applied to series circuits.

<p align="center"><img  alt="figure 3" src="https://electrica.ual.es/spacor/images/spacorcasethree.png" width=50% height=50%></p>

In this circuit it is satisfied that,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=v = Ri + L \frac{di}{dt}"></p>
 
In contrast to the parallel model, it is now convenient to define the geometric vector as,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=z = i\sigma_1 + i'\sigma_2 + v\sigma_3"></p>

so that we work with the information contained in the derivative of the current instead of the integral of the voltage. Developing the vector z and taking into account KVL, this yield,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=z = i \underbrace{(\sigma_1 + R \sigma_3)}_{e_1} + i'\underbrace{(\sigma_2 + L\sigma_3)}_{e_2}"></p>
  
The plane that expands <img src="https://render.githubusercontent.com/render/math?math=e_1"> and <img src="https://render.githubusercontent.com/render/math?math=e_3"> is,
Therefore the relationship between the bivectors  and  gives  and the relationship between  and  gives . Let us now proceed with the calculation of the oscillating plane,
The  and  values are,
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=K = e_1 \wedge e_2 = (\sigma_1 + R\sigma_3)\wedge(\sigma_2 + L\sigma_3)"></p>
 
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=K = \sigma_{12} - R\sigma_{23} + L\sigma_{13}"></p>
 
[Matlab simulation for case 3](examples/matlabCase3.md)
