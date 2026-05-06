## Case 4: Mixed series and parallel circuit with sinusoidal source

The technique can also be used to solve mixed series and parallel circuits. For example, be the case of a resistor <img src="https://render.githubusercontent.com/render/math?math=R_p">, in parallel with a real air-core coil of values <img src="https://render.githubusercontent.com/render/math?math=R_s"> y <img src="https://render.githubusercontent.com/render/math?math=L">  as in the figure. Kirchhoff's equations for this circuit tell us that

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=i = i_R + i_L \rightarrow i_L = i - i_R = i -Gv"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=v = R_S i_L + L i"></p>

Working on the two previous equations, we can obtain the value of the total current as a function of the voltage, its derivative and the derivative of the current,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=v = R_s(i-Gv) + L(i' - Gv')"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=i = \frac{1 + R_sG}{R_s}v- \frac{L}{R_s}i' + \frac{LG}{R_s}v'"></p>

<p align="center"><img  alt="figure 3" src="https://electrica.ual.es/spacor/images/spacorcasefour.png" width=50% height=50%></p> 

The new vector has the following components,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=y = v\sigma_1 + v'\sigma_2 +i'\sigma_3 + i \sigma_4"></p>

For this case we add the information of the derivative of the current and voltage so that we can compute the 3 unknown parameters. Developing the above expression we arrive at,
 
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=y = v\underbrace{(\sigma_1+\frac{1+R_sG}{R_s}\sigma_4)}_{e_1}+v'\underbrace{(\sigma_2+\frac{LG}{R_s}\sigma_4)}_{e_2}+i'\underbrace{(\sigma-\frac{L}{R_s}\sigma_4)}_{e_3}"></p>
 
We calculate the 3D subspace,
 
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=K = e_1 \wedge e_2 \wedge e_3"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=K = \sigma_{123} - \frac{L}{R_s}\sigma_{124} - \frac{LG}{R_s}\sigma{134}+\frac{1+R_sG}{R_s}\sigma_{234} "></p> 
 
In this case, <img src="https://render.githubusercontent.com/render/math?math=G"> is obtained as the ratio between <img src="https://render.githubusercontent.com/render/math?math=\sigma_{134}"> and <img src="https://render.githubusercontent.com/render/math?math=\sigma_{124}">. Having found <img src="https://render.githubusercontent.com/render/math?math=G">, we perform the quotient between <img src="https://render.githubusercontent.com/render/math?math=\sigma_{234}">  and <img src="https://render.githubusercontent.com/render/math?math=\sigma_{123}"> (which we will call a), so that we obtain 
<img src="https://render.githubusercontent.com/render/math?math=R_s = \frac{1}{a-G}">. Having found <img src="https://render.githubusercontent.com/render/math?math=R_s"> , we make the ratio between <img src="https://render.githubusercontent.com/render/math?math=\sigma_{124}"> and <img src="https://render.githubusercontent.com/render/math?math=\sigma_{123}"> and multiply by **Rs** and obtain **L**. Finally, we compute the oscillating trivector with the help of the tangent and normal vectors,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=y = v\sigma_1 + v'\sigma_2 + i'\sigma_3 + i \sigma_4"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=y' = v'\sigma_1 + v''\sigma_2 + i''\sigma_3 + i' \sigma_4"></p> 
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=y ''= v''\sigma_1 + v'''\sigma_2 + i'''\sigma_3 + i'' \sigma_4"></p>

In this way, the 3D volume (osculant) is computed as,

<p align="center"><img src="https://render.githubusercontent.com/render/math?math=K_{osc} = y \wedge y' \wedge y''"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=K_{osc} = [i'''(vv''-v'^2)-v'''(vi''-i'v')+v''(v'i''-i'v'')]\sigma_{123}"></p> 
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=+ [i''(vv''-v'^2)-v'''(vi'-iv')+v''(v'i'-iv'')]\sigma_{124}"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=+[i''(vi''-i'v')-i'''(vi'-iv')+v''(i'^2-ii'')]\sigma_{134} "></p> 
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=+[i''(v'i''-i'v'')-i'''(v'i'-iv'')+v'''(i'^2-ii'')]\sigma_{234} "></p>

Thus, the values of the elements are as follows,
 
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=G = \frac{i''(vi''-i'v')-i'''(vi'-iv')+v''(i'^2-ii'')}{i''(vv''-v'^2)-v'''(vi'-iv')+v''(v'i'-iv'')} "></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=a = \frac{i''(vi''-i'v')-i'''(vi'-iv')+v''(i'^2-ii'')}{i'''(vv''-v'^2)-v'''(vi''-i'v')+v''(v'i''-i'v'')}"></p> 
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=R_s = \frac{1}{a-G}"></p>
<p align="center"><img src="https://render.githubusercontent.com/render/math?math=L = -R_s\frac{i''(vv''-v'^2)-v'''(vi'-iv')+v''(v'i'-iv'')}{i'''(vv''-v'^2)-v'''(vi''-i'v')+v''(v'i''-i'v'')}"></p> 


[Matlab simulation for case 4](examples/matlabCase4.md)
