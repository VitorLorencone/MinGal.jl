# Installation

```@setup set
using Mingal
Algebra(3)
```

*How to install*

This package is not yet in [Metadata.jl](https://github.com/JuliaLang/METADATA.jl), so the package can't be installed with the Julia package manager. Thus, for installing it, from the Julia REPL, type ] to enter the Pkg REPL mode and run:

```@julia
pkg> add https://github.com/VitorLorencone/Mingal.jl.git
```

## First Steps

Now, for start Mingal, return to Julia mode in REPL and type:

```julia
julia> using Mingal
```

The second step to use Mingal is define an Algebra, an environment. The environment is defined through the `Algebra()` function such as:

```julia
julia> Algebra(3, 0)
```

In this case, were created the 3D space. More information about the created space, is just showed in the REPL. By now, try to execute the following command:

```julia
julia> (1+e1e2e3)*(e1)
1.0*e1 + 1.0*e2e3
```