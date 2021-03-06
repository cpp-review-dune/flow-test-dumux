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
 * \brief Exercise to show the diffusive spreading of contaminants.
 */

#ifndef DUMUX_LENS_2P2C_PROPERTIES_HH
#define DUMUX_LENS_2P2C_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>

#include <dumux/material/fluidsystems/h2on2.hh>

#include "lens2pspatialparams.hh"
#include "lens2p2cproblem.hh"

//////////
// Specify the properties for the lens problem
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct LensTwoPTwoCProblem { using InheritsFrom = std::tuple<TwoPTwoC, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::LensTwoPTwoCProblem> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::LensTwoPTwoCProblem> { using type = LensTwoPTwoCProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LensTwoPTwoCProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LensTwoPTwoCProblem>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Lens2pSpatialParams<FVGridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
