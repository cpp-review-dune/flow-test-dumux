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
 * \brief application, using the 3p box model
 */
#include "config.h"
#include "properties.hh"

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>



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
                    errorMessageOut += "\n\nThe list of mandatory options for this program is:\n"
                                        "\t-TimeManager.TEnd              End of the simulation [s] \n"
                                        "\t-TimeManager.DtInitial         Initial timestep size [s] \n"
                                        "\t-Grid.UpperRight               Coordinates of upper right grid corner [m] \n"
                                        "\t-Grid.Cells                    Number of cells in x, y direction [-] \n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv) try
{

    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::InfiltrationThreePBoxTypeTag;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();


    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(fvGridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))  // Fehlt ebenso
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);


    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<FVGridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    //! set some check points for the time loop
    timeLoop->setPeriodicCheckPoint(getParam<Scalar>("TimeLoop.EpisodeLength"));

    // time loop
    timeLoop->start(); do
    {
        // set the end of the next time step as time in the problem to control
        // the boundary conditions for the implicit Euler scheme
        problem->setTime(timeLoop->time()+timeLoop->timeStepSize());

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // write vtk output on check points
        if (timeLoop->isCheckPoint() || timeLoop->timeStepIndex() == 1 || timeLoop->finished())
            vtkWriter.write(timeLoop->time());

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main

catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                "). Most likely, the DGF file name is wrong "
                "or the DGF file is corrupted, "
                "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
