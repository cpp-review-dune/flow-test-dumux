// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/


#ifndef DUMUX_LENS2P_PROPERTIES_HH
#define DUMUX_LENS2P_PROPERTIES_HH


#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/material/components/tabulatedcomponent.hh>

#include "lens2pspatialparams.hh"
#include "lens2pproblem.hh"



//////////
// Specify the properties for the lens problem
//////////

namespace Dumux::Properties {


// Create new type tags
namespace TTag {
struct LensTwoPProblem { using InheritsFrom = std::tuple<TwoP, BoxModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::LensTwoPProblem> { using type = LensTwoPProblem<TypeTag>; };


// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::LensTwoPProblem> { using type = Dune::YaspGrid<2>; };



// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LensTwoPProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::TabulatedComponent<Components::H2O<Scalar>> >;
    using NonwettingPhase = FluidSystems::OnePGas<Scalar, Components::N2<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LensTwoPProblem>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Lens2pSpatialParams<FVGridGeometry, Scalar>;
};


} // end namespace Dumux::Properties

#endif
