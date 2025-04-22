# Examples

```@setup ss
using Mingal
Complex = Algebra(0, 1, ["i"])
```

*How to use?*

!!! note
    Since the system is in its implementation phase, new examples will be added in the future.

## Complex Plane

This example demonstrates a way to work with complex numbers using the Mingal library, in a user-friendly manner with regard to notation:

```julia
Complex = Algebra(0, 1, ["i"])
```

```@repl ss
i*i == -1*id
5*i*i*i == -5*i
(1+i)*(1+i) == 2*i
(1+i)*(1-i) == 2*id
```

```@setup xx
using Mingal
R2 = Algebra(2,0,["i","j"])
```

## ℝ² plane

This example demonstrates a way to work with the ℝ² vectors using the Mingal library, in a user-friendly manner with regard to notation:

```julia
R2 = Algebra(2,0,["i","j"])
```

```@repl xx
i^i == 0*id
i\j == 0*id
i^j == ij
(i + 2*j)\(5*i+2*j) == 9*id
```