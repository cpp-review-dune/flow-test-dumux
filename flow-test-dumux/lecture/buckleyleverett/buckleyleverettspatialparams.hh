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
#ifndef BUCKLEYLEVERETT_SPATIALPARAMS_HH
#define BUCKLEYLEVERETT_SPATIALPARAMS_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

//forward declaration
template<class TypeTag>
class BuckleyLeverettSpatialParams;

namespace Properties {

// The spatial parameters TypeTag
namespace TTag {
struct BuckleyLeverettSpatialParamsTypeTag {};
}

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BuckleyLeverettSpatialParamsTypeTag> { using type = BuckleyLeverettSpatialParams<TypeTag>; };

} // end namespace Properties

// forward declaration
class LinearMaterialDefault;
class LinearMaterial;

template<class TypeTag>
class BuckleyLeverettSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using CoordScalar =  typename Grid::ctype;
    static constexpr int dim=Grid::dimension;
    static constexpr int dimWorld=Grid::dimensionworld, numEq=1;
    using Element = typename Grid::Traits::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using FieldMatrix = Dune::FieldMatrix<Scalar,dim,dim>;
    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:

    BuckleyLeverettSpatialParams(const Problem& problem)
    : ParentType(problem)
    , pcKrSwCurve_("SpatialParams")
    {
        Scalar permFactor = 1.0;  //0.001/(1000*9.81);
        constPermeability_ = getParam<Scalar>("SpatialParams.Permeability")*permFactor;
        porosity_ = getParam<Scalar>("SpatialParams.Porosity");
    }

    static constexpr bool pcSwCurveIsLinear()
    { return (std::is_same_v<PcKrSwCurve, LinearMaterial> || std::is_same_v<PcKrSwCurve, LinearMaterialDefault>); }

    Scalar intrinsicPermeability (const Element& element) const
    { return constPermeability_; }

    Scalar porosity(const Element &element) const
    { return porosity_; }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return pcKrSwCurve_; }

private:
    const PcKrSwCurve pcKrSwCurve_;
    Scalar constPermeability_;
    Scalar porosity_;
    Scalar swr_;
    Scalar snr_;
    Scalar pcEntry_;
    Scalar lambda_;

};

} // end namespace Dumux

#endif
