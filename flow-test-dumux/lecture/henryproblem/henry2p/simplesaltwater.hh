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
 * \brief A much simpler (and thus potentially less buggy) version of
 *        pure water.
 */
#ifndef DUMUX_SIMPLE_SALTWATER_HH
#define DUMUX_SIMPLE_SALTWATER_HH

#include <cmath>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/components/componenttraits.hh>

#include <dumux/material/idealgas.hh>

namespace Dumux {
/*!
 *
 * \brief A much simpler (and thus potentially less buggy) version of
 *        pure water.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class SimpleSaltwater
: public Components::Base<Scalar, SimpleSaltwater<Scalar> >
, public Components::Liquid<Scalar, SimpleSaltwater<Scalar> >
, public Components::Gas<Scalar, SimpleSaltwater<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;

    static const Scalar R;  // specific gas constant of water

public:
    /*!
     * \brief A human readable name for the water.
     */
    static std::string name()
    {
        return "H2O";
    }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of water.
     */
    static Scalar molarMass()
    {
        return 18e-3;
    }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of water.
     */
    static Scalar criticalTemperature()
    {
        return 647.096; /* [K] */
    }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of water.
     */
    static Scalar criticalPressure()
    {
        return 22.064e6; /* [N/m^2] */
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at water's triple point.
     */
    static Scalar tripleTemperature()
    {
        return 273.16; /* [K] */
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at water's triple point.
     */
    static Scalar triplePressure()
    {
        return 611.657; /* [N/m^2] */
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar vaporPressure(Scalar T)
    {
        if (T > criticalTemperature())
            return criticalPressure();

        if (T < tripleTemperature())
            return 0; // water is solid: We don't take sublimation into account

        static const Scalar n[10] =
        {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        Scalar sigma = T + n[8]/(T - n[9]);
        Scalar A = (sigma + n[0])*sigma + n[1];
        Scalar B = (n[2]*sigma + n[3])*sigma + n[4];
        Scalar C = (n[5]*sigma + n[6])*sigma + n[7];

        Scalar tmp = Scalar(2.0)*C/(std::sqrt(B*B - Scalar(4.0)*A*C) - B);
        tmp *= tmp;
        tmp *= tmp;

        return Scalar(1e6)*tmp;
    }

    /*!
     * \brief Specific enthalpy of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        return 1976*(temperature - 293.15) + 2.45e6;
    }

    /*!
     * \brief Specific enthalpy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    {
        return 4180*(temperature - 293.15);
    }

    /*!
     * \brief Specific internal energy of steam \f$\mathrm{[J/kg]}\f$.
     *
     *        Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     *
     *        Rearranging for internal energy yields: \f$u = h - pv\f$.
     *
     *        Exploiting the Ideal Gas assumption (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     *        The universal gas constant can only be used in the case of molar formulations.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            1/molarMass()* // conversion from [J/(mol K)] to [J/(kg K)]
            IdealGas::R*temperature; // = pressure *spec. volume for an ideal gas
    }

    /*!
     * \brief Specific internal energy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    {
        return
            liquidEnthalpy(temperature, pressure) -
            pressure/liquidDensity(temperature, pressure);
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    {
        return true;
    }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    {
        return false;
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of steam at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return molarMass()*IdealGas::molarDensity(temperature, pressure);
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    {
        return true;
    }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief The density of pure water at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1025.;
    }

    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    {
        return liquidDensity(temperature, pressure)/molarMass();
    }

    /*!
     * \brief The pressure of water in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The liquid pressure is undefined for incompressible fluids");
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure, bool regularize = true)
    {
        return 1e-05;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 1e-03;
    }

    /*!
     * \brief Specific isobaric heat capacity of the component [J/(kg K)] as a liquid.
     *
     *        \param temperature temperature of component in \f$\mathrm{[K]}\f$
     *        \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *        source: http://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&Action=Page
     *        @ T= 281.15 (8°C) , p=0.1MPa)
     */
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        return 4.2e3 ;
    }
};

template <class Scalar>
const Scalar SimpleSaltwater<Scalar>::R = Dumux::Constants<Scalar>::R / 18e-3;

} // end namespace Dumux

#endif
