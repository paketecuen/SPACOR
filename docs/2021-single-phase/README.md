# SPACOR Space Vectors for Parameter Identification of Non-Linear Circuits

It proposes a systematic method for the identification of the load circuit parameters (say the R, L, and C elements) based only on the information of the instantaneous volt age and current measured at the point of common coupling (pcc). **Geometric Algebra (GA)** and concepts of differential geometry are used to produce a rigorous mathematical framework. The identification is formulated as a multidimensional geometrical problem that is solved conveniently by means of GA. Once the passive elements of the load have been identified, the active and reactive powers can be computed from first electromagnetic principles (Maxwell Equations). The theory is general and is verified with linear and nonlinear circuits. The paper shows single-phase circuits but the theory can be extended to three phase circuits. The method is easy to program and has shown to be very robust for all tested cases. Because of its generality, the method presented will find applications beyond electric circuits.
https://electrica.ual.es/spacor

A discussion about it can be found in the paper "Determination of Instantaneous Powers from a Novel Time-Domain Parameter Identification Method of Non-Linear Single-Phase Circuits. Francisco G. Montoya, Francisco De Leon,  Franscico M. Arrabal-Campos, and Alfredo Alcayde. IEEE Transactions on Power Delivery. 2021"
https://doi.org/10.1109/TPWRD.2021.3133069

```bibtext
@ARTICLE{9640464,
  author={Montoya, Francisco G and De Leon, Francisco and Arrabal-Campos, Francisco M. and Alcayde, Alfredo},
  journal={IEEE Transactions on Power Delivery}, 
  title={Determination of Instantaneous Powers from a Novel Time-Domain Parameter Identification Method of Non-Linear Single-Phase Circuits}, 
  year={2021},
  volume={},
  number={},
  pages={1-1},
  doi={10.1109/TPWRD.2021.3133069}}
``` 

## EXAMPLES

[Case 1: Parallel RL load with Sinusoidal source](Case1.md)

[Case 2: Parallel RL load switching the resitive part with Sinusoidal source](Case2.md)

[Case 3: Serie RL with Sinusoidal source](Case3.md)

[Case 4: Mixed series and parallel circuit with sinusoidal source](Case4.md)

[Case 5: Parallel RL with piecewise non-linear inductor](Case5.md)

[Case 6: Parallel RLC with switching capacitor and resistor](Case6.md)

[Case 7: Series RLC with switching capacitor and resistor](Case7.md)

 
