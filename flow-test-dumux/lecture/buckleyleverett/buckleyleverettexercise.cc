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
 * \brief DOC ME!
 */
#include "config.h"
#include "buckleyleverettproperties.hh"
#include <dumux/common/start.hh>
#include <dumux/common/properties.hh>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe List of Mandatory arguments for this program is:\n"
                                        "\t-TimeManager.TEnd                                The end of the simulation. [s] \n"
                                        "\t-TimeManager.DtInit                              The initial timestep size. [s] \n"
                                        "\t-SpatialParams.Permeability                      The intrinsic permeability of the porous medium [m^2]\n"
                                        "\t-SpatialParams.Porosity                          The porosity of the porous medium [-]\n"
                                        "\t-SpatialParams.BrooksCoreyLambda                 The pore size distribution parameter for the \n"
                                        "\t                                                 \t Brooks-Corey capillary pressure - saturation relationship [-]\n"
                                        "\t-SpatialParams.BrooksCoreyPcEntry                The entry pressure for the \n"
                                        "\t                                                 \t Brooks-Corey capillary pressure - saturation relationship [Pa]\n"
                                        "\t-SpatialParams.Swr                               The residual saturation of the wetting phase [-]\n"
                                        "\t-SpatialParams.Snr                               The residual saturation of the nonwetting phase [-]\n"
                                        "\t-Fluid.DensityW                                  The density of the wetting phase [kg/m^3]\n"
                                        "\t-Fluid.DensityNW                                 The density of the nonwetting phase [kg/m^3]\n"
                                        "\t-Fluid.ViscosityW                                The dynamic viscosity of the wetting phase [kg/(ms)]\n"
                                        "\t-Fluid.ViscosityNW                               The dynamic viscosity of the nonwetting phase [kg/(ms)]\n"
                                        "\t-Grid.NumberOfCellsX                             The grid resolution in x direction [-]\n"
                                        "\n optional: \n"
                                        "\t-Output.ParaviewOutput                           Boolean, default is not writing ViPLab but paraview output";
        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
    using TypeTag = Dumux::Properties::TTag::BuckleyLeverettProblemTypeTag;
    return Dumux::start<TypeTag>(argc, argv, usage);
}
