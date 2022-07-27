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
/**
 * \file
 * \brief properties file for Henry exercise
  *
 */
#ifndef DUMUX_HENRY1P2C_PROPERTIES_HH
#define DUMUX_HENRY1P2C_PROPERTIES_HH

#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>

#include <dune/grid/yaspgrid.hh>

#include "henry1p2cspatialparameters.hh"
#include "henry1p2cproblem.hh"

#include "watersaltfluidsystem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Henry1p2cProblemTypeTag { using InheritsFrom = std::tuple<BoxModel, OnePNC>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Henry1p2cProblemTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Henry1p2cProblemTypeTag> { using type = Henry1p2cProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Henry1p2cProblemTypeTag> { using type = WaterSaltFluidSystem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Henry1p2cProblemTypeTag>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Henry1p2cSpatialParams<FVGridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableCompositionalDispersion<TypeTag, TTag::Henry1p2cProblemTypeTag> { static constexpr bool value = true; };

template<class TypeTag>
struct CompositionalDispersionModel<TypeTag, TTag::Henry1p2cProblemTypeTag> { using type = ScheideggersDispersionTensor<TypeTag>; };

//Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Henry1p2cProblemTypeTag> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
