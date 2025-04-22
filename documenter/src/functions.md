# Functions

```@setup ss
using Mingal
```

*What are the main functions and what they do*

## Algebra()

Algebra() is the most important function of any system with this library, as it defines an environment, that is, an algebra that will be used to handle future operations. Thus, as its definition follows, we have:

```@docs
Mingal.Algebra
```

Which allows us to create spaces like the following, representing the imaginary set:

```@example ss
Algebra(0,1,["i"])
```

## Blades

All types introduced in Mingal are subtype of `AbstractGeometricAlgebraType`. Blade type is the most primitive subtype wich is used to define more complete structures, like Multivectors. Once we define an Algebra, a set of Blade is created. For example, let us create the space R3.

```@example ss
Algebra(3)
```

Now, the elements `id`, `e1`, `e2`, `e3`, `e1e2`, `e1e3`, `e2e3`, `e1e2e3` were created, representing the basis blades of the Algebra.

```@example ss
typeof(e1) <: Mingal.AbstractGeometricAlgebraType
```

By now, operation tables related to the geometric, internal and outer products are created and will be used in these operations. The `id` element is a BasisBlade used in place of the number 1.

## Multivectors

Multivectors are mutable structures primarily used for internal operations within Mingal, but they are subtypes of `AbstractGeometricAlgebraType`. They are composed of the sum of blades arranged in a sparse vector, used for more efficient representation.

Note that there when is defined a scalar product or sum or diference between an real scalar and BasisBlade element a new Multivector is always returned. A Multivector struct is showed just in function and order of the Basis Blade set that define the space.

```@repl ss; continued=true
Multivector <: Mingal.AbstractGeometricAlgebraType # true
Multivectors([8],[2.5]) == 2.5*e1e2e3 # true
Multivectors([1, 2, 3, 4, 5, 6],[1, 2, 3, 4, 5, 6]) == 1 + 2*e1 + 3*e2 + 4*e3 + 5*e1e2 + 6*e1e3 # true
```

## Simple Blade functions

### bladeindex()

```@docs
Mingal.bladeindex
```

And it works as the following example:

```@repl ss
Mingal.bladeindex(e1e2)
```

### bladescalar()

```@docs
Mingal.bladescalar
```

And it works as the following example:

```@repl ss
Mingal.bladescalar(4.5*e1e2)
```

### grade()

```@docs
Mingal.grade
```

And it works as the following example:

```@repl ss
grade(e1e2e3)
```

### gradeprojection()

```@docs
Mingal.gradeprojection
```

And it works as the following example:

```@repl ss
gradeprojection(-e1e2, 2) == -1.0*e1e2
```

## Operations

### Scalar Multiplication

The scalar multiplication is a type of operation that occurs between a blade or multivector and any scalar number. Its symbol for operations is `*`

And it works as the following example:

```@repl ss
2*id*(e1+e2) == 2.0*e1+2.0*e2
```

### Geometric Product

The geometric product is the core of all geometric algebra, and here in Mingal, it is represented between blades and multivectors using the symbol `*`, which should not be confused with the scalar product symbol. Although they are the same, they have different contexts and operands of different types.

And it works as the following example:

```@repl ss
e3*e1 == -e1e3
```

### Inner Product

The inner product is another type of operation defined for Geometric Algebra and is represented within Mingal by the symbol `\`.

And it works as the following example:

```@repl ss
e1e2\e1 == -e2
```

### Outer Product

The Outer product is another type of operation defined for Geometric Algebra and is represented within Mingal by the symbol `^`.

And it works as the following example:

```@repl ss
e1^e1e2 == 0*id
```

### All Together

The products can be used in whichever way is preferred, including all mixed together, with sums, parentheses, or any usual arithmetic operation.

And it works as the following example:

```@repl ss
(1+2*e1)\(5e1e2^e2) == 0.0*id
```