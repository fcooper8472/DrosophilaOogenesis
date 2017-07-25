/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#ifndef TESTOOGENESIS_HPP_
#define TESTOOGENESIS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellAncestorWriter.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

// Header files included in this project
#include "GermariumBoundaryCondition.hpp"

// Should usually be called last.
#include "PetscSetupAndFinalize.hpp"

class TestOogenesis : public AbstractCellBasedTestSuite
{
public:

    void Test3dOogenesis() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        /* First we set up some numbers that will define the crypt and villus geometry */
        double germarium_radius = 1.0;

        /* We then create a couple of cells at the base of each germarium.
         * (we put two cells in each crypt to set off delta-notch patterning) */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  0.0, 0.0, 0.0));

        /*
         * We then convert this list of nodes to a `NodesOnlyMesh`,
         * which doesn't do very much apart from keep track of the nodes.
         */
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        /*
         * Next we have to create the cells that will be associated with these nodes.
         * So we make an empty vector in which to store the cells and then loop over
         * each node, adding cells as we go.
         */
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        CellsGenerator<OocyteAggregateCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements());

        /*
         * We now create a cell population, which keeps track of a mesh and cells and the association between them.
         * In this case we need a `NodeBasedCellPopulation` in three dimensions.
         */
        NodeBasedCellPopulation<3> germarium(mesh, cells);

        /* We limit the absolute movement that cells can make to cause error messages if numerics become unstable */
        germarium.SetAbsoluteMovementThreshold(10);

        /* We then instruct the cell population to output some useful information for plotting in VTK format in e.g. paraview */
        germarium.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        /*
         * We now set up our cell-based simulation class.
         * We have sub-classed the main simulator class `OffLatticeSimulation`
         * which is in the core of Chaste with a new simulation class `SimplifiedDeltaNotchOffLatticeSimulation`
         * which can be found in this project's `src` folder.
         *
         * The reason for this is that this class performs calculations of average Delta levels surrounding each cell
         * at the end of each timestep. You will see a minimum number of methods have been overridden, and the class
         * is fairly simple.
         */
        OffLatticeSimulation<3> simulator(germarium);
        simulator.SetOutputDirectory("DrosophilaOogenesis");
        simulator.SetDt(1.0/120.0);
        /* We limit the output to every 120 time steps (1 hour) to reduce output file sizes */
        simulator.SetSamplingTimestepMultiple(120);

        /*
         * We now create a force law and pass it to the simulation
         * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
         */
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(30.0); // default is 15.0;
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        /*
         * The most complex part of this problem definition is that of the boundary condition that limits
         * cell locations to a 2D surface in 3D space. This has been defined in a separate class
         * `MultipleCryptGeometryBoundaryCondition` which can be found in this project's `src` folder.
         */
        auto p_boundary_condition = boost::make_shared<GermariumBoundaryCondition>(args);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /* We then set an end time and run the simulation */
        simulator.SetEndTime(250.0);
        simulator.Solve(); // to 250 hours
    }
};

#endif /*TESTOOGENESIS_HPP_*/
