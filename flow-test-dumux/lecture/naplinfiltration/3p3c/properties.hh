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
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */
#ifndef DUMUX_NAPLINFILTRATIONPROPERTIES_3P_3C_HH
#define DUMUX_NAPLINFILTRATIONPROPERTIES_3P_3C_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3p3c/model.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>

#include "myh2oairmesitylene.hh"
#include "../spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties
{

// Create new type tags
namespace TTag {
struct InfiltrationThreePThreeCTypeTag { using InheritsFrom = std::tuple<ThreePThreeC>; };
struct InfiltrationThreePThreeCBoxTypeTag { using InheritsFrom = std::tuple<InfiltrationThreePThreeCTypeTag, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::InfiltrationThreePThreeCTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::InfiltrationThreePThreeCTypeTag> { using type = InfiltrationThreePThreeCProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::InfiltrationThreePThreeCTypeTag>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InfiltrationSpatialParams<FVGridGeometry, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::InfiltrationThreePThreeCTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::MyH2OAirMesitylene<Scalar>;
};

} // end namespace Dumux::Properties

#endif
