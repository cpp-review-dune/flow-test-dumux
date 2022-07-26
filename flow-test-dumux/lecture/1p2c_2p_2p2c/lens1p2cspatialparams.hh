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
 * \brief Class for defining spatial parameters
 */
#ifndef DUMUX_1P2C_SPATIALPARAMS_HH
#define DUMUX_1P2C_SPATIALPARAMS_HH

#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

//forward declaration
template<class FVGridGeometry, class Scalar>
class Lens1p2cSpatialParams;

namespace Properties {
// The spatial parameters TypeTag
namespace TTag {
struct Lens1p2cSpatialParams {};
}

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Lens1p2cSpatialParams>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Lens1p2cSpatialParams<FVGridGeometry, Scalar>;
};

} // end namespace Properties

/*!
 * \brief Class for defining spatial parameters
 */
template<class FVGridGeometry, class Scalar>
class Lens1p2cSpatialParams : public FVPorousMediumFlowSpatialParamsOneP<FVGridGeometry, Scalar,  Lens1p2cSpatialParams<FVGridGeometry, Scalar> >
{
    using ThisType = Lens1p2cSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     */
    Lens1p2cSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        lensLowerLeft_[0] = 0.0;
        lensLowerLeft_[1] = 0.0;
        lensUpperRight_[0]= 4.1;
        lensUpperRight_[1]= 1.0;

        lensPorosity_ = getParam<double>("SpatialParams.Fine.Porosity");
        outerPorosity_ = getParam<double>("SpatialParams.Coarse.Porosity");

        lensK_ = getParam<double>("SpatialParams.Fine.Permeability");
        outerK_ = getParam<double>("SpatialParams.Coarse.Permeability");
        temperature_ = getParam<double>("SpatialParams.Temperature");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return permeability
     */
    template<class ElementSolution>
    Scalar permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        if (isInLens_(scv.dofPosition()))
            return lensK_;

        else
            return outerK_;
    }

    /*!
     * \brief Function for defining the porosity.
     *        That is possibly solution dependent.
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        if (isInLens_(scv.dofPosition()))
            return lensPorosity_;

        else
            return outerPorosity_;
    }

    /*!
     * \brief set the lens coordinates
     */
    void setLensCoords(const GlobalPosition& lensLowerLeft, const GlobalPosition& lensUpperRight)
    {
        lensLowerLeft_ = lensLowerLeft;
        lensUpperRight_ = lensUpperRight;
    }

    /*!
     * \brief Returns the temperature at the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return temperature_;
    }

private:
    bool isInLens_(const GlobalPosition &pos) const
    {
        for (int i = 0; i < dimWorld; ++i)
        {
            if (pos[i] < lensLowerLeft_[i] || pos[i] > lensUpperRight_[i])
                return false;
        }

        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;
    Scalar lensPorosity_;
    Scalar outerPorosity_;
    Scalar temperature_;
};

} // end namespace Dumux

#endif
