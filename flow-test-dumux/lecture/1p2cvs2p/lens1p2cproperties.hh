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
 * \brief Soil decontamination problem where a contaminant infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENS_1P2C_PROPERTIES_HH
#define DUMUX_LENS_1P2C_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/discretization/box.hh>
#include "../1p2c_2p_2p2c/lens1p2cspatialparams.hh"
#include "lens1p2cproblem.hh"


template <typename TypeTag>
class LensOnePTwoCProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct LensOnePTwoCProblem { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::LensOnePTwoCProblem> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::LensOnePTwoCProblem> { using type = LensOnePTwoCProblem<TypeTag>; };

// set fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LensOnePTwoCProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/false>>;
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LensOnePTwoCProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
public:
    using type = Lens1p2cSpatialParams<FVGridGeometry, Scalar>;
};

} //end namespace Dumux::Properties

#endif
