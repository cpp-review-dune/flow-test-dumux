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
#ifndef DUMUX_LECTURE_MM_BUCKLEYLEVERETT_BUCKLEYLEVERETTPROPERTIES_HH
#define DUMUX_LECTURE_MM_BUCKLEYLEVERETT_BUCKLEYLEVERETTPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>

//pseudo-oil and h2o have to be used to make viscosity and density input parameters
#include "pseudooil.hh"
#include "pseudoh2o.hh"
#include "buckleyleverettspatialparams.hh"

//the problem file
#include "buckleyleverettproblem.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct BuckleyLeverettProblemTypeTag { using InheritsFrom = std::tuple<BuckleyLeverettSpatialParamsTypeTag, IMPESTwoP, FVTransportTwoP, FVPressureTwoP>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::BuckleyLeverettProblemTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::BuckleyLeverettProblemTypeTag> { using type = BuckleyLeverettProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BuckleyLeverettProblemTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, PseudoOil<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, PseudoH2O<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

template<class TypeTag>
struct EvalCflFluxFunction<TypeTag, TTag::BuckleyLeverettProblemTypeTag> { using type = EvalCflFluxCoats<TypeTag>; };
template<class TypeTag>
struct Formulation<TypeTag, TTag::BuckleyLeverettProblemTypeTag> { static constexpr int value = SequentialTwoPCommonIndices::pwsw; };

} // end namespace Dumux::Properties

#endif
