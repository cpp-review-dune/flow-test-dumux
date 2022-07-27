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
/*!
 * \file
 *
 * \brief Properties file of soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_HENRY2PPROPERTIES_HH
#define DUMUX_HENRY2PPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include "simplesaltwater.hh"

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/2p/model.hh>

#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "henry2pspatialparams.hh"
#include "henry2pproblem.hh"

//////////
// Specify the properties for the Henry2p problem
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Henry2pProblemTypeTag { using InheritsFrom = std::tuple<BoxModel, TwoP>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Henry2pProblemTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Henry2pProblemTypeTag> { using type = Henry2pProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Henry2pProblemTypeTag>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Henry2pSpatialParams<FVGridGeometry, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Henry2pProblemTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, SimpleSaltwater<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

} // end namespace Dumux::Properties

#endif
