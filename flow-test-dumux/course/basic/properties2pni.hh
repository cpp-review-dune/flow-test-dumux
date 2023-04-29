// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief The two-phase nonisothermal porousmediumflow properties file for exercise-basic
 */

#ifndef DUMUX_EX_BASIC_PROPERTIES_2PNI_HH
#define DUMUX_EX_BASIC_PROPERTIES_2PNI_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "injection2pniproblem.hh"
#include "injection2pspatialparams.hh"

namespace Dumux::Properties {

 /*!
* TODO:dumux-course-task 4
* Inherit from the TwoPNI model instead of TwoP here
*/
// Create new type tags
namespace TTag {
struct Injection2pNITypeTag { using InheritsFrom = std::tuple<TwoP>; };
struct Injection2pNICC { using InheritsFrom = std::tuple<Injection2pNITypeTag, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Injection2pNITypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Injection2pNITypeTag> { using type = Injection2PNIProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Injection2pNITypeTag>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = InjectionSpatialParams<GridGeometry, Scalar>;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Injection2pNITypeTag>
{
    using type = FluidSystems::H2ON2< GetPropType<TypeTag, Properties::Scalar>,
                                      FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true> >;
};

} // end namespace Dumux::Properties

#endif
