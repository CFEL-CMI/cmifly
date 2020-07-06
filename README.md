# CMI fly

A simple Python program for trajectory-based simulations of the deflection of molecular beams in
inhomogeneous electric fields.

The initial version of this code was published in the following paper and any use of the code should
reference that publication:

> Yuan-Pin Chang, Daniel Horke, Sebastian Trippel, Jochen Küpper: Spatially-controlled complex
> molecules and their applications, [_International Reviews in Physical Chemistry_ **34**, 557–590
> (2015)](https://dx.doi.org/10.1080/0144235x.2015.1077838); [arXiv:1505.05632
> [physics]](https://arxiv.org/abs/1505.05632).

See the accompanying [license](./LICENSE.md) and the [documentation](#documentation) for further
details.


## Getting Started with CMIfly

For now, get the latest version from [github](https://github.com/CFEL-CMI/cmifly) and run
`CMIfly.py`.


## Documentation

For now, this code is undocumented; this needs to be changed.

### Some notes on the usage for different setups

For the a-type deflector, use the following lines in CMIfly.py:
```
deflector_fieldnorm_filename     = 'deflector_field_norm'
deflector_fieldgradient_filename = 'deflector_field_gradient'
# deflector_voltage for which the fields were calculated (kV)
deflector_field_voltage = 5.
# geometric boundary of the deflector
rod_center = [0, 0.0036] # rod center (x,y)
rod_radius = 0.003 # rod radius
trough_center = [0, 0.002408] # trough center (x,y)
trough_radius = 0.0032 # trough radius
```

For the b-type deflector, use the following lines in CMIfly.py:
```
deflector_fieldnorm_filename     = 'Ec.norm.txt'
deflector_fieldgradient_filename = 'gt.grad.txt'
# deflector_voltage for which the fields were calculated (kV)
deflector_field_voltage = 60.
# geometric boundary of the deflector
rod_center = [0, 0.0036] # rod center (x,y)
rod_radius = 0 #0.003 # rod radius
trough_center = [0, 0.002408] # trough center (x,y)
trough_radius = 0 #0.0032 # trough radius
Since the rod_radius = 0, it will now use the fields to determine the edge of the deflector
```



<!-- Put Emacs local variables into HTML comment
Local Variables:
coding: utf-8
fill-column: 100
End:
-->
