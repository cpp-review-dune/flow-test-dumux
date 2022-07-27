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
#ifndef DUMUX_BUCKLEYLEVERETT_ANALYTICAL_HH
#define DUMUX_BUCKLEYLEVERETT_ANALYTICAL_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>
#include <dumux/material/fluidmatrixinteractions/2p/materiallaw.hh>

/*!
 * \file
 * \brief  Analytical solution of the buckley-leverett problem
 * \author Markus Wolff
 */

namespace Dumux {
/*!
 * \brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * the Buckley-Leverett problem
 */

template<class TypeTag> class BuckleyLeverettAnalytic
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using CellData = GetPropType<TypeTag, Properties::CellData>;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimworld = GridView::dimensionworld;
    static constexpr int wPhaseIdx = Indices::wPhaseIdx;
    static constexpr int nPhaseIdx = Indices::nPhaseIdx;
    static constexpr int satEqIdx = Indices::satEqIdx;
    using BlockVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1> >;
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    // functions needed for analytical solution

    void initializeAnalytic()
    {
        analyticSolution_.resize(size_);
        analyticSolution_ = 0;
        error_.resize(size_);
        error_ = 0;
        elementVolume_.resize(size_);
        elementVolume_ = 0;
    }

    /*!
     * \brief DOC ME!
     * \param materialLawParams DOC ME!
     */
    void prepareAnalytic()
    {
        const auto& dummyElement = *problem_.gridView().template begin<0>();
        const auto fluidMatrixInteraction = problem_.spatialParams().fluidMatrixInteractionAtPos(dummyElement.geometry().center());
        swr_ = fluidMatrixInteraction.effToAbsParams().swr();
        snr_ = fluidMatrixInteraction.effToAbsParams().snr();
        porosity_ = problem_.spatialParams().porosity(dummyElement);
        time_ = 0;
        satVec_ = swr_;

        for (int i = 1; i < pointNum_; i++)
        {
            satVec_[i] = satVec_[i - 1] + (1 - snr_ - swr_) / intervalNum_;
        }

        FluidState fluidState;
        fluidState.setTemperature(problem_.temperatureAtPos(dummyGlobal_));
        fluidState.setPressure(wPhaseIdx, problem_.referencePressureAtPos(dummyGlobal_));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressureAtPos(dummyGlobal_));
        Scalar viscosityW = FluidSystem::viscosity(fluidState, wPhaseIdx);
        Scalar viscosityNW = FluidSystem::viscosity(fluidState, nPhaseIdx);

        for (int i = 0; i < pointNum_; i++)
        {
            fractionalW_[i] = fluidMatrixInteraction.krw(satVec_[i])/viscosityW;
            fractionalW_[i] /= (fractionalW_[i] + fluidMatrixInteraction.krn(satVec_[i])/viscosityNW);
        }

        dfwdsw_ = 0;
        for (int i = 1; i < intervalNum_; i++)
        {
            dfwdsw_[i] = (fractionalW_[i + 1] - fractionalW_[i - 1]) / (satVec_[i + 1] - satVec_[i - 1]);
        }

        for (int i = 0; i < pointNum_; i++)
        {
            if (dfwdsw_[i] > dfwdsw_[i + 1])
            {
                dfwdswmax_ = i;
                break;
            }
        }
    }

    void calcSatError()
    {
        error_ = 0;
        elementVolume_ = 0;
        ElementIterator eItEnd = problem_.gridView().template end<0> ();

        for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
        {
            // get entity
            const Element& element = *eIt;
            int index = problem_.variables().index(*eIt);
            elementVolume_[index] = element.geometry().volume();
        }

        for (int i = 0; i < size_; i++)
        {
            error_[i] = analyticSolution_[i] - problem_.variables().cellData(i).saturation(wPhaseIdx);
        }
    }

    void updateExSol()
    {
        // position of the fluid front
        xf_ = 0;
        for (int i = 0; i < pointNum_; i++)
        {
            xf_[i] = vTot_ * time_ / porosity_ * dfwdsw_[i];
        }

        // position of maximum xf_
        int xfmax = 0;
        int xhelp = pointNum_ / 3;
        int xhelpold = 0;
        int xhelpoldold = 0;
        int xhelp2 = 0;
        for (int i = 0; i < pointNum_; i++)
        {
            if (xf_[i] > xf_[i + 1])
            {
                xfmax = i;
                break;
            }
        }

        // balancing of the areas ahead of the front and below the curve
        bool a = true;
        Scalar A1;
        Scalar A2;
        Scalar b;

        while (a)
        {
            A1 = 0;

            for (int i = 0; i <= xhelp - 1; i++)
            {
                A1 += (satVec_[i] - swr_ + satVec_[i + 1] - swr_) * 0.5 * (xf_[i + 1] - xf_[i]);
            }

            A2 = 0;

            for (int i = xhelp; i <= xfmax - 1; i++)
            {
                A2 += (satVec_[xfmax] - satVec_[i] + satVec_[xfmax] - satVec_[i + 1]) * 0.5 * (xf_[i + 1] - xf_[i]);
            }

            b = xf_[xfmax];
            xhelp2 = xfmax;

            while (b > xf_[xhelp])
            {
                xhelp2 += 1;
                b = xf_[xhelp2];
            }

            for (int i = xfmax; i <= xhelp2; i++)
            {
                A2 += (satVec_[i] - satVec_[xfmax] + satVec_[i + 1] - satVec_[xfmax]) * 0.5 * (xf_[i] - xf_[i + 1]);
            }

            xhelpoldold = xhelpold;
            xhelpold = xhelp;

            if (fabs(A1) > fabs(A2))
            {
                xhelp = xhelp - 1;
            }

            if (fabs(A1) < fabs(A2))
            {
                xhelp = xhelp + 1;
            }

            if (xhelp == xhelpoldold)
            {
                a = false;
            }
        }

        // iterate over vertices and get analytic saturation solution
        ElementIterator eItEnd = problem_.gridView().template end<0> ();
        for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
        {
            // get global coordinate of cell center
            const GlobalPosition& globalPos = eIt->geometry().center();

            // find index of current vertex
            int index = problem_.variables().index(*eIt);

            // account for linear material law
            if constexpr (SpatialParams::pcSwCurveIsLinear())
            {
                if (globalPos[0] <= xf_[1])
                {
                    analyticSolution_[index] = 1 - snr_;
                }
                if (globalPos[0] > xf_[1])
                {
                    analyticSolution_[index] = swr_;
                }
            }
            // non-linear material law
            else
            {
                // find x_f next to global coordinate of the vertex
                int xnext = 0;
                for (int i = intervalNum_; i >= 0; i--)
                {
                    if (globalPos[0] < xf_[i])
                    {
                        xnext = i;
                        break;
                    }
                }

                // account for the area not yet reached by the front
                if (globalPos[0] > xf_[xhelp2])
                {
                    analyticSolution_[index] = swr_;
                    continue;
                }

                if (globalPos[0] <= xf_[xhelp2])
                {
                    analyticSolution_[index] = satVec_[xnext];
                    continue;
                }
            }
        }

        // call function to calculate the saturation error
        calcSatError();
    }

    void calculateAnalyticSolution()
    {
        time_ += problem_.timeManager().timeStepSize();
        updateExSol();
    }

    BlockVector AnalyticSolution() const
    { return analyticSolution_; }

    //Write saturation and pressure into file
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        BlockVector *analyticSolution = writer.allocateManagedBuffer (size_);
        BlockVector *error = writer.allocateManagedBuffer (size_);

        *analyticSolution = analyticSolution_;
        *error = error_;

        writer.attachCellData(*analyticSolution, "saturation (exact solution)");
        writer.attachCellData(*error, "error_");
    }

    //! Construct an IMPES object.
    BuckleyLeverettAnalytic(Problem& problem, Scalar totalVelocity = 3e-7)
        : problem_(problem)
        , analyticSolution_(0)
        , error_(0)
        , elementVolume_(0)
        , size_(problem.gridView().size(0))
        , vTot_(totalVelocity)
    {
        dummyGlobal_ = 0.0;
        dummyGlobal_[0] = 1.0;
    }

    void initialize()
    {
        initializeAnalytic();
        prepareAnalytic();
    }

private:
    Problem& problem_;

    BlockVector analyticSolution_;
    BlockVector error_;
    BlockVector elementVolume_;

    Scalar time_;
    int size_;
    Scalar swr_;
    Scalar snr_;
    Scalar porosity_;
    Scalar vTot_;
    static constexpr int intervalNum_ = 1000, pointNum_ = intervalNum_ + 1;
    Dune::FieldVector<Scalar, pointNum_> satVec_;
    Dune::FieldVector<Scalar, pointNum_> fractionalW_;
    Dune::FieldVector<Scalar, pointNum_> dfwdsw_;
    Dune::FieldVector<Scalar, pointNum_> xf_;
    int dfwdswmax_;
    GlobalPosition dummyGlobal_;

};

} // end namespace Dumux

#endif
