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

#include <properties/proliferative_types/TransitCellProliferativeType.hpp>
#include "GermariumBoundaryCondition.hpp"

#include "CellLabel.hpp"
#include "DrosophilaOogenesisEnumerations.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"

#include "Debug.hpp"

GermariumBoundaryCondition::GermariumBoundaryCondition(AbstractCellPopulation<3>* pCellPopulation,
                                                       double germariumRadius)
        : AbstractCellPopulationBoundaryCondition<3>(pCellPopulation),
          mRadiusOfGermarium(germariumRadius)
{
    if (dynamic_cast<NodeBasedCellPopulation<3>*>(this->mpCellPopulation) == NULL)
    {
        EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
    }
}

double GermariumBoundaryCondition::GetRadiusOfGermarium() const
{
    return mRadiusOfGermarium;
}

void GermariumBoundaryCondition::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    // Iterate over the cell population
    for (auto cell_iter = this->mpCellPopulation->Begin(); cell_iter != this->mpCellPopulation->End(); ++cell_iter)
    {
        const c_vector<double,3> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

        // We get the cell label from the current cell and obtain its colour
        auto abstract_property = cell_iter->rGetCellPropertyCollection().GetProperties<CellLabel>().GetProperty();
        unsigned cell_colour = boost::static_pointer_cast<CellLabel>(abstract_property)->GetColour();

        // If the cell is an aggregate, it stays on the x axis
        if (cell_colour == TYPE_AGGREGATE)
        {
            // We map the y and z locations to zero and just keep the x location
            c_vector<double,3> new_location = zero_vector<double>(3);

            // The stem cell will remain at the origin; other nodes are just mapped down to the x axis
            if(not cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
            {
                new_location[0] = cell_location[0];
            }

            // Move node to new location
            p_node->rGetModifiableLocation() = new_location;
        }
        // If it's a follicle cell, it is more free to move
        else if (cell_colour == TYPE_FOLLICLE)
        {
            // The stem cells will remain at their starting point; other nodes are just mapped down to the x axis
            if(cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
            {
                p_node->rGetModifiableLocation()[0] = 0.0;
                p_node->rGetModifiableLocation()[2] = 0.0;

                p_node->rGetModifiableLocation()[1] = p_node->rGetLocation()[1] > 0.0 ? 1.0 : -1.0;
            }
            else
            {
                // Ensure node cannot be in negative x
                p_node->rGetModifiableLocation()[0] = cell_location[0] < 0.0 ? 0.0 : cell_location[0];

                c_vector<double, 2> y_z_location;
                y_z_location[0] = cell_location[1];
                y_z_location[1] = cell_location[2];

                // If the node is outside the radius of the germarium, map it back to the surface (keeping x the same)
                if (norm_2(y_z_location) > mRadiusOfGermarium)
                {
                    y_z_location /= norm_2(y_z_location);

                    p_node->rGetModifiableLocation()[1] = y_z_location[0];
                    p_node->rGetModifiableLocation()[2] = y_z_location[1];
                }
            }
        }
        else
        {

        }
    }
}

bool GermariumBoundaryCondition::VerifyBoundaryCondition()
{
    return true;
}


void GermariumBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RadiusOfGermarium>" << mRadiusOfGermarium << "</RadiusOfGermarium>\n";
    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<3>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(GermariumBoundaryCondition)
