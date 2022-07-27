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
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 */
#ifndef DUMUX_PLUMESHAPE_PROBLEM_HH
#define DUMUX_PLUMESHAPE_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/discretization/box/scvftoscvboundarytypes.hh>

#include <lecture/common/co2tablesbenchmark3.hh>

// per default use isothermal model
#ifndef ISOTHERMAL
#define ISOTHERMAL 0
#endif

namespace Dumux {

/*!
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 *
 * The domain is sized 200m times 100m and consists of four layers, a
 * permeable reservoir layer at the bottom, a barrier rock layer with reduced permeability, another reservoir layer
 * and at the top a barrier rock layer with a very low permeablility.
 *
 * CO2 is injected at the permeable bottom layer
 * from the left side. The domain is initially filled with brine.
 *
 * The grid is unstructered and permeability and porosity for the elements are read in from the grid file. The grid file
 * also contains so-called boundary ids which can be used assigned during the grid creation in order to differentiate
 * between different parts of the boundary.
 * These boundary ids can be imported into the problem where the boundary conditions can then be assigned accordingly.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is false.
 *
 * To run the simulation execute the following line in shell (works with the box and cell centered spatial discretization method):
 * <tt>./test_ccco2 </tt> or <tt>./test_boxco2 </tt>
 */
template <class TypeTag >
class PlumeShapeProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    // copy some indices for convenience
    enum
    {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        // phase presence index
        firstPhaseOnly = Indices::firstPhaseOnly,

        // component indices
        BrineIdx = FluidSystem::BrineIdx,
        CO2Idx = FluidSystem::CO2Idx,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,
        contiCO2EqIdx = conti0EqIdx + CO2Idx
    };

#if !ISOTHERMAL
    enum {
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
    };
#endif

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using CO2 = Components::CO2<Scalar, GeneratedCO2Tables::CO2Tables>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

    // the discretization method we are using
    static constexpr auto discMethod = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;

    // world dimension to access gravity vector
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    PlumeShapeProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        nTemperature_         = getParam<int>("FluidSystem.NTemperature");
        nPressure_            = getParam<int>("FluidSystem.NPressure");
        pressureLow_          = getParam<Scalar>("FluidSystem.PressureLow");
        pressureHigh_         = getParam<Scalar>("FluidSystem.PressureHigh");
        temperatureLow_       = getParam<Scalar>("FluidSystem.TemperatureLow");
        temperatureHigh_      = getParam<Scalar>("FluidSystem.TemperatureHigh");
        name_                 = getParam<std::string>("Problem.Name");
        pressure_             = getParam<Scalar>("InitialConditions.Pressure"); // hydrodynamic pressure at top layer
        temperature_          = getParam<Scalar>("InitialConditions.Temperature");
        injectionRate_        = getParam<Scalar>("BoundaryConditions.InjectionRate");
        injectionTemperature_ = getParam<Scalar>("BoundaryConditions.InjectionTemperature");
        dipAngle_             = getParam<Scalar>("SpatialParams.DipAngle"); // [deg]

        dipAngleRadians_ = (dipAngle_ * M_PI)/180.0; // [rad]
        // rotate the coordinate system / gravity:
        gravity_[1] = -9.81 * cos(dipAngleRadians_);
        gravity_[0] = -9.81 * sin(dipAngleRadians_);


        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        //stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position
     *
     * This problem assumes a geothermal gradient with
     * a surface temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     */

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllNeumann();

        // set East boundaries to Dirichlet values
        if ( globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_  )
        {
            values.setAllDirichlet();
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param returns the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub control volume face
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector fluxes(0.0);
        // kg/(m^2*s) or mole/(m^2*s) depending on useMoles
        const auto globalPos = scvf.ipGlobal();
        if ( globalPos[1] < 30.0 + eps_ && globalPos[0] < eps_ )
        {
            Scalar massFlux = injectionRate_ /30.0; // [kg/(s*m^2)]
            fluxes[contiCO2EqIdx] = -massFlux; // kg/(s*m^2)
#if !ISOTHERMAL
            const Scalar pressure = elemVolVars[scvf.insideScvIdx()].pressure(FluidSystem::phase1Idx);
            fluxes[energyEqIdx] = -massFlux*CO2::gasEnthalpy(injectionTemperature_, pressure);
#endif
        }

        return fluxes;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values at a position
     *
     * \returns the initial values for the conservation equations in
     *           \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * The internal method for the initial condition
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(firstPhaseOnly);

        const Scalar densityW = FluidSystem::Brine::liquidDensity(temperature_, 1e5);

        const Scalar moleFracLiquidCO2 = 0.00;
        const Scalar moleFracLiquidBrine = 1.0 - moleFracLiquidCO2;

        const Scalar meanM = FluidSystem::molarMass(BrineIdx)*moleFracLiquidBrine +
                             FluidSystem::molarMass(CO2Idx)*moleFracLiquidCO2;

        const Scalar massFracLiquidCO2 = moleFracLiquidCO2*FluidSystem::molarMass(CO2Idx)/meanM;

        // set initial pressure and pressure Dirichlet conditions:
        const Scalar distanceToTopBoundary_ = ( this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1] ) - globalPos[1]; // distance without rotation
        const Scalar distanceToTopBoundaryRotated_ = distanceToTopBoundary_ / std::cos(dipAngleRadians_);

        const Scalar distanceTopBoundaryToSurface_ = std::sin(dipAngleRadians_) * (this->gridGeometry().bBoxMax()[0] - globalPos[0]);
        values[Indices::pressureIdx] = pressure_ + densityW * 9.81 * distanceToTopBoundaryRotated_ +
                                       densityW * 9.81 * distanceTopBoundaryToSurface_;

        values[Indices::switchIdx] = massFracLiquidCO2;

#if !ISOTHERMAL
        values[temperatureIdx] = temperature_;//290.0 + (depthBOR_ - globalPos[1])*0.03;
#endif


        return values;
    }

    Scalar injectionRate_;

#if !ISOTHERMAL
    Scalar injectionTemperature_;
#endif

    int nTemperature_;
    int nPressure_;

    Scalar temperature_;
    Scalar pressure_;

    Scalar dipAngleRadians_, dipAngle_;

    GlobalPosition gravity_;

    Scalar eps_ = 1.e-6;

    std::string name_ ;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
};

} //end namespace Dumux

#endif
