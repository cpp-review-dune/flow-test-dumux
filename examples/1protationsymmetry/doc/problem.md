<!-- Important: This file has been automatically generated by generate_example_docs.py. Do not edit this file directly! -->


| [:arrow_left: Back to the main documentation](../README.md) | [:arrow_right: Continue with part 2](main.md) |
|---|---:|

# Part 1: Simulation setup

The code documentation is structured as follows:

[[_TOC_]]


## Compile-time settings (`properties.hh`)
This file defines the `TypeTag` used for the simulation in this example, for
which we specialize a number of compile-time `properties`.

<details open>
<summary><b>Click to hide/show the file documentation</b> (or inspect the [source code](../properties.hh))</summary>

### Includes
<details><summary> Click to show includes</summary>

```cpp

#include <dune/grid/yaspgrid.hh> // for `Dune::YaspGrid`
#include <dumux/discretization/box.hh> // for `TTag::BoxModel`
```

The `OneP` type tag specializes most of the `properties` required for single-
phase flow simulations in DuMu<sup>x</sup>. We will use this in the following to inherit the
respective properties, and subsequently specialize those properties for our
type tag, which we want to modify or for which no meaningful default can be set.

```cpp
#include <dumux/porousmediumflow/1p/model.hh> // for `TTag::OneP`
```

The local residual for incompressible flow is included.
The one-phase flow model (included above) uses a default implementation of the
local residual for single-phase flow. However, in this example we are using an
incompressible fluid phase. Therefore, we are including the specialized local
residual which contains functionality to analytically compute the entries of
the Jacobian matrix. We will use this in the main file.

```cpp
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>
```

We will use a single liquid phase consisting of a component with constant fluid properties.

```cpp
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
```

As mentioned at the beginning of the documentation of this example, DuMu<sup>x</sup>
provides specialized implementations of the extrusion type for rotational symmetric geometries.
The extrusion class takes care of adjusting volume and area
computations to account for the extrusion about the symmetry axes.

```cpp
#include <dumux/discretization/extrusion.hh>
```

The classes that define the problem and parameters used in this simulation

```cpp
#include "problem.hh"
#include "spatialparams.hh"
```

</details>

### `TypeTag` definition
A `TypeTag` for our simulation is defined, which inherits properties from the
single-phase flow model and the box scheme.

```cpp
namespace Dumux::Properties {
namespace TTag {
struct OnePRotSym { using InheritsFrom = std::tuple<OneP, BoxModel>; };
}
```

### Property specializations

In the following piece of code, mandatory `properties` for which no meaningful
default can be set, are specialized for our type tag `OnePRotSym`.

```cpp
// We use a structured 1D grid with an offset. This allows us to define the
// computational domain to be between the radii $`r_1`$ and $`r_2`$ as illustrated
// in the beginning of the documentation of this example
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePRotSym>
{ using type =  Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>; };

// The problem class specifying initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePRotSym>
{ using type = RotSymExampleProblem<TypeTag>; };

// Our spatial parameters class defining the permeability and porosity of the porous medium:
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePRotSym>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = RotSymExampleSpatialParams<GridGeometry, Scalar>;
};

// We use a single liquid phase consisting of a component with constant fluid properties.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePRotSym>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};
```

As mentioned before, we modify the areas and volumes used
for integrals over the control volume and the control volume faces by changing the extrusion type.
Here, we pass these traits to the grid geometry of the box scheme (the scheme
that we use here) and specialize the `GridGeometry` property accordingly.

```cpp
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::OnePRotSym>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    // We take the default traits as basis and exchange the extrusion type
    // The first axis (x-axis) is the radial axis, hence the zero. That means we rotate about the second axis (y-axis).
    struct GGTraits : public BoxDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

public:
    // Pass the above traits to the box grid geometry such that it uses the
    // rotation-symmetric sub-control volumes and faces.
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};
```

Moreover, here we use a local residual specialized for incompressible flow
that contains functionality related to analytic differentiation.

```cpp
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePRotSym>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties
```


</details>



## The problem class (`problem.hh`)
This file contains the __problem class__ which defines the initial and boundary
conditions for the single-phase flow simulation.

<details open>
<summary><b>Click to hide/show the file documentation</b> (or inspect the [source code](../problem.hh))</summary>

### Includes

```cpp
#include <cmath> // for `std::log`
#include <dumux/common/boundarytypes.hh> // for `BoundaryTypes`
#include <dumux/common/properties.hh> // for `GetPropType`
#include <dumux/common/parameters.hh> // for `getParam`
#include <dumux/porousmediumflow/problem.hh>  // for `PorousMediumFlowProblem`
```

### The problem class
We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.

```cpp
namespace Dumux {

template<class TypeTag>
class RotSymExampleProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
```

In the constructor, we obtain a number of parameters, related to fluid
properties and boundary conditions, from the input file.

```cpp
    RotSymExampleProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // fluid properties
        k_ = getParam<Scalar>("SpatialParams.Permeability");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        rho_ = getParam<Scalar>("Component.LiquidDensity");

        // The inner radius r1 can be determined from the grid
        r1_ = gridGeometry->bBoxMin()[0];

        // boundary conditions
        q1_ = getParam<Scalar>("Problem.Q1"); // mass flux into the domain at r1 in kg/s/m
        p1_ = getParam<Scalar>("Problem.P1"); // pressure at the inner boundary at r1
    }
```

#### Specify the types of boundary conditions
This function is used to define the type of boundary conditions used depending on the location.
Two types of boundary  conditions can be specified: Dirichlet or Neumann boundary condition.
On a Dirichlet boundary, the values of the primary variables need to be fixed. On a Neumann
boundary condition, values for derivatives need to be fixed. Here, we use Dirichlet boundary
conditions on all boundaries.

```cpp
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }
```

#### Specify Dirichlet boundary condition values
This function is used to specify the values of the primary variables at Dirichlet boundaries.
Here, we evaluate the analytical solution (see below) to define the pressures at the boundaries.

```cpp
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return exactSolution(globalPos); }
```

#### Analytical solution
The analytical solution to the problem of this example reads:

```math
p = p (r) = p_1 - \frac{q_1 \nu}{2 \pi k} \text{ln} (\frac{r}{r_1}),
```

where $`q_1`$ is the mass flux into the domain at the inner radius $`r_1`$
(in kg/s/m) and $`\nu = \mu/\varrho`$ is the kinematic viscosity.
The following function evaluates this solution depending on the
position in the domain. We use this function here both to specify Dirichlet
boundaries and to evaluate the error of the numerical solutions obtained for
different levels of grid refinement.

```cpp
    PrimaryVariables exactSolution(const GlobalPosition& globalPos) const
    {
        const auto r = globalPos[0];
        const auto p = p1_ - 1.0/(2*M_PI)*nu_/k_*q1_*std::log(r/r1_);
        return p;
    }

    const Scalar exactVelocity(const GlobalPosition& globalPos) const
    {
        const auto r = globalPos[0];
        const auto v = q1_/(2*M_PI)/rho_/r;
        return v;
    }

private:
    // private data members required for the analytical solution
    Scalar q1_, k_, nu_, r1_, p1_, rho_;
};

} // end namespace Dumux
```


</details>



## Parameter distributions (`spatialparams.hh`)

This file contains the __spatial parameters class__ which defines the
distributions for the porous medium parameters permeability and porosity
over the computational grid.

<details open>
<summary><b>Click to hide/show the file documentation</b> (or inspect the [source code](../spatialparams.hh))</summary>

We include the spatial parameters class for single-phase models discretized
by finite volume schemes, from which the spatial parameters defined for this
example inherit.

```cpp
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
```

### The spatial parameters class

In the `RotSymExampleSpatialParams` class, we define the functions needed to describe
the porous medium, that is, porosity and permeability.
We inherit from the `FVPorousMediumFlowSpatialParamsOneP` class here, which is the base class
for spatial parameters in the context of single-phase porous medium flow
applications using finite volume discretization schemes.

```cpp
namespace Dumux {

template<class GridGeometry, class Scalar>
class RotSymExampleSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, RotSymExampleSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RotSymExampleSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    // Spatial parameter classes for porous medium flow applications need to
    // export the type used for intrinsic permeabilities.
    using PermeabilityType = Scalar;

    // In the constructor we obtain the permeability value from the input file.
    RotSymExampleSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { permeability_ = getParam<Scalar>("SpatialParams.Permeability"); }
```

#### Porosity distribution
This function is used to define the porosity distribution in the
computational domain. Here, we use a constant porosity of 1.0.

```cpp
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }
```

#### Permeability distribution
This function is used to define the permeability distribution in the
computational domain. Here, we use a constant permeability that is
defined in the input file.

```cpp
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

private:
    Scalar permeability_;
};

} // end namespace Dumux
```


</details>


| [:arrow_left: Back to the main documentation](../README.md) | [:arrow_right: Continue with part 2](main.md) |
|---|---:|

