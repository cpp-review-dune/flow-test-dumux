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
 * \brief The one-phase porousmediumflow problem for exercise mainfile
 */
#ifndef DUMUX_EX_MAINFILE_PROPERTIES_HH
#define DUMUX_EX_MAINFILE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>


#include <dumux/porousmediumflow/1p/model.hh>
// TODO: dumux-course-task 3
// uncomment the incompressiblelocalresidual which is a specialization of the standard immisible localresidual for one phase incompressible cases and provides an analytic jacobian.
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "1pspatialparams.hh"
#include "1pproblem.hh"

namespace Dumux::Properties {

// Create the new type tag nodes:
// Here we define the incompressible type tag as well as the compressible type tag.
// The incompressible uses a different fluidsystem than the compressible
namespace TTag {
struct OnePBase { using InheritsFrom = std::tuple<OneP>; };
struct OnePIncompressible { using InheritsFrom = std::tuple<OnePBase, CCTpfaModel>; };
struct OnePCompressible { using InheritsFrom = std::tuple<OnePBase, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePBase> { using type = Dune::YaspGrid<2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePBase> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePBase> {
    using type = OnePTestSpatialParams<GetPropType<TypeTag, GridGeometry>, GetPropType<TypeTag, Scalar>>;
};

// the fluid system for incompressible tests
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePIncompressible>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// TODO: dumux-course-task 3
// set the OneP Incompressible local residual for the OnePIncompressible type tag. This provides an analytic jacobian to be used for the analytic solution. Change that by setting:
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePIncompressible> { using type = OnePIncompressibleLocalResidual<TypeTag>; };


// the fluid system for compressible tests
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePCompressible>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::TabulatedComponent<Components::H2O<Scalar>>>;
};

// Disable caching (for testing purposes)
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePBase> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePBase> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePBase> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
