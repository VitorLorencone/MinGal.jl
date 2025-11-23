<p align="center">
  <img src="docs/assets/logo-dark.png" alt="Logo" width="50%">
</p>

# MinGal

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://vitorlorencone.github.io/MinGal.jl/)

*A Geometric Algebra library written in Julia*

This page is dedicated to document the use of [MinGal](https://github.com/VitorLorencone/MinGal.jl), a library for [Geometric Algebra](https://en.wikipedia.org/wiki/Geometric_algebra). Our goal is to provide an accessible and easy to learn experience for students and researchers interested in using (and programming with) a Geometric Algebra library.

Our aim is to offer a relatively simple, open, easy-to-install, and focused environment for Julia users. This library was developed by implementing basic operations as described in textbooks, without an emphasis on optimizing functions based on other libraries. For the future, the plans are to expand the ease of use and introduce concepts for modeling and visualizing spaces in Geometric Algebra.

This project is continually evolving, and many aspects are being refined to achieve better results.

Check the [documentation](https://vitorlorencone.github.io/MinGal.jl/) for more information on how to use, examples and other possibilities.

## Developed By

This project is being developed as part of a scientific initiation program at the State University of Maringá.

- [Vitor Madeira Lorençone](https://github.com/VitorLorencone) - Student
- [Emerson Vitor Castelani](https://github.com/evcastelani) - Professor

## Instalattion

This package is not yet registered, so it can't be installed with the Julia package manager. Thus, for installing it, from the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add https://github.com/VitorLorencone/MinGal.jl
```

Or you could also run:

```julia
julia> Using Pkg
julia> Pkg.add(url="https://github.com/VitorLorencone/MinGal.jl")
```

## First Steps

Now, to start MinGal, return to Julia mode in REPL and type:

```julia
julia> using MinGal
```

The second step to use Mingal is define an Algebra, an environment. The environment is defined through the `Algebra()` function such as:

```julia
julia> Algebra(3)
```

In this case, we created the 3D space. More information about the created space, is just showed in the REPL. By now, try to execute the following command:

```julia
julia> (1+e1e2e3)*(e1)
1.0*e1 + 1.0*e2e3
```

Don't forget to check the [documentation](https://vitorlorencone.github.io/MinGal.jl/) for more information.