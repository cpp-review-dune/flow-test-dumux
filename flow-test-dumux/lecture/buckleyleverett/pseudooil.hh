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
 * \brief Properties of pure water \f$H_2O\f$.
 */
#ifndef DUMUX_PSEUDOOIL_HH
#define DUMUX_PSEUDOOIL_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>

namespace Dumux {

/*!
 * \brief Rough estimate for testing purposes of some oil.
 */
template <class ScalarT>
class PseudoOil : public Components::Base<ScalarT, PseudoOil<ScalarT> >
                , public Components::Liquid<ScalarT, PseudoOil<ScalarT> >
{
public:
    using Scalar = ScalarT;

    /*!
     * \brief A human readable name for the water.
     */
    static std::string name()
    {
        return "Oil";
    }

    /*!
     * \brief Rough estimate of the density of oil [kg/m^3].
     * \param temperature DOC ME!
     * \param pressure DOC ME!
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return density_;
    }

    /*!
     * \brief Rough estimate of the viscosity of oil kg/(ms).
     * \param temperature DOC ME!
     * \param pressure DOC ME!
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return viscosity_;
    }

    /*!
     * \brief DOC ME!
     * \param viscosity DOC ME!
     */
    static void setViscosity(Scalar viscosity)
    {
        viscosity_ =  viscosity;
    }

    /*!
     * \brief DOC ME!
     * \param density DOC ME!
     */
    static void setDensity(Scalar density)
    {
        density_ =  density;
    }

private:
    static Scalar viscosity_;
    static Scalar density_;
};

template <class ScalarT>
typename PseudoOil<ScalarT>::Scalar PseudoOil<ScalarT>::viscosity_ = 0.01;

template <class ScalarT>
typename PseudoOil<ScalarT>::Scalar PseudoOil<ScalarT>::density_ = 1000;


} // end namepace Dumux

#endif
