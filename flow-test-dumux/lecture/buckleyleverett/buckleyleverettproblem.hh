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
#ifndef DUMUX_LECTURE_MM_BUCKLEYLEVERETT_BUCKLEYLEVERETTPROBLEM_HH
#define DUMUX_LECTURE_MM_BUCKLEYLEVERETT_BUCKLEYLEVERETTPROBLEM_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>

#include "buckleyleverettanalyticsolution.hh"
#include <lecture/common/viplaboutput.hh>

//pseudo-oil and h2o have to be used to make viscosity and density input parameters
#include "pseudooil.hh"
#include "pseudoh2o.hh"

namespace Dumux {

template<typename TypeTag>
class BuckleyLeverettProblem: public IMPESProblem2P<TypeTag>
{
    using ParentType = IMPESProblem2P<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using PrimaryVariables = typename GetProp<TypeTag, Properties::SolutionTypes>::PrimaryVariables;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int wPhaseIdx = Indices::wPhaseIdx;
    static constexpr int nPhaseIdx = Indices::nPhaseIdx;
    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int swIdx = Indices::swIdx;
    static constexpr int pressEqIdx = Indices::pressureEqIdx;
    static constexpr int satEqIdx = Indices::satEqIdx;
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief DOC ME!
     * \param timeManager DOC ME!
     * \param gridView DOC ME!
     */
    BuckleyLeverettProblem(TimeManager& timeManager, typename GridView::Grid &grid)
    : ParentType(timeManager, grid)
    , eps_(1e-8)
    , pLeftBc_(2.0e5)
    , analyticSolution_(*this, 3e-7)
    , viplabOutput_(*this)
    {
        // check the cell size restriction
        const auto cells = getParam<std::vector<int> >("Grid.Cells");
        if (cells[1] != 1)
            DUNE_THROW(Dune::IOError, "This is a 1D problem. Do not change the number of cells in y-direction!");

        if (cells[0] > 200)
            DUNE_THROW(Dune::IOError, "Maximum allowed number of cells in x-direction for this problem is: 200!");

        upperRight_ = getParam<GlobalPosition>("Grid.UpperRight");

        PseudoOil<Scalar>::setViscosity( getParam<Scalar>("Fluid.ViscosityW") );
        PseudoH2O<Scalar>::setViscosity( getParam<Scalar>("Fluid.ViscosityNW") );

        PseudoOil<Scalar>::setDensity( getParam<Scalar>("Fluid.DensityW") );
        PseudoH2O<Scalar>::setDensity( getParam<Scalar>("Fluid.DensityNW") );

        densityNonwetting_ = getParam<Scalar>("Fluid.DensityNW");

        swr_ = getParam<Scalar>("SpatialParams.Swr");
        snr_ = getParam<Scalar>("SpatialParams.Snr");

        paraviewOutput_ = getParam<bool>("Output.paraviewOutput", true);
    }

    void init()
    {
        ParentType::init();
        analyticSolution_.initialize();
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
    const std::string name() const
    {
        return "buckleyleverett";
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    void postTimeStep()
    {
        analyticSolution_.calculateAnalyticSolution();
        ParentType::postTimeStep();
    }

    void addOutputVtkFields()
    {
        ParentType::addOutputVtkFields();
        analyticSolution_.addOutputVtkFields(this->resultWriter());
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     * \param globalPos DOC ME!
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 273.15 + 10; // -> 10°C
    }

    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1.0e5; // -> 1 bar
    }
    /*!
     * \brief DOC ME!
     * \param values DOC ME!
     * \param globalPos DOC ME!
     */
    void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
    {
        values = 0;
    }
    /*!
     * \brief DOC ME!
     * \param bcTypes DOC ME!
     * \param globalPos DOC ME!
     */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
    {
        if (globalPos[0] < eps_)//all Dirichlet at west boundary
        {
            bcTypes.setAllDirichlet();
        }

        else if (globalPos[0] > upperRight_[0] - eps_)//at east boundary
        {
            bcTypes.setNeumann(pressEqIdx);//Neumann for the pressure equation
            bcTypes.setOutflow(satEqIdx);//Outflow for the transport equation
        }

        // all other boundaries, i.e., top and bottom, Neumann
        else
        {
            bcTypes.setAllNeumann();
        }
    }
    /*!
     * \brief DOC ME!
     * \param values DOC ME!
     * \param globalPos DOC ME
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values[pressureIdx] = pLeftBc_;
        values[swIdx] = 1.0 - snr_;
    }
    /*!
     * \brief DOC ME!
     * \param values DOC ME!
     * \param globalPos DOC ME
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0; //no-flow Neumann at top and bottom

        if (globalPos[0] > upperRight_[0] - eps_) //east boundary
        {
            // the volume flux should remain constant, when density is changed
            // here, we multiply by the density of the Nonwetting Phase
            const Scalar referenceDensity = 1000.0;
            values[nPhaseIdx] = 3e-4 * densityNonwetting_/referenceDensity;
        }
    }
    /*!
     * \brief DOC ME!
     * \param values DOC ME!
     * \param globalPos DOC ME
     */
    void initialAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = pLeftBc_;
        values[swIdx] = swr_;
    }

    /**
     * \brief Override outputfunction for ViPLab output
     */
    void writeOutput()
    {
        if (paraviewOutput_)
            ParentType::writeOutput(false);

        else
        {
            std::size_t numElements = this->gridView().size(0);
            viplabOutput_.writeHeader();
            std::vector<int> color({0, 0, 255}); // blue
            std::vector<Scalar> solutionVectorSaturation(numElements);
            for (std::size_t i = 0; i < numElements; i++)
                solutionVectorSaturation[i] = this->variables().cellData(i).saturation(wPhaseIdx);

            viplabOutput_.writeVipLabOutput(solutionVectorSaturation, color);

            for (std::size_t i = 0; i < numElements; i++)
                solutionVectorSaturation[i]=analyticSolution_.AnalyticSolution()[i];

            color = {255, 0, 0}; // red
            viplabOutput_.writeVipLabOutput(solutionVectorSaturation, color);
            viplabOutput_.writeLegend();
        }
    }

private:
    GlobalPosition upperRight_;
    Scalar eps_, swr_, snr_;
    Scalar pLeftBc_;
    Scalar densityNonwetting_;
    BuckleyLeverettAnalytic<TypeTag> analyticSolution_;
    ViplabOutput<TypeTag> viplabOutput_;
    bool paraviewOutput_;
};

} // end namespace Dumux

#endif
