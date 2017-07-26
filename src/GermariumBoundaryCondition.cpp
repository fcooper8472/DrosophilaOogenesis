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

#include "GermariumBoundaryCondition.hpp"

#include "DrosophilaOogenesisEnumerations.hpp"
#include "NodeBasedCellPopulation.hpp"

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

        // We map the y and z locations to zero and just keep the x location
        c_vector<double,3> new_location = zero_vector<double>(3);

        // Node with ID zero will remain at the origin; other nodes are just mapped down to the x axis
        if(cell_iter->GetCellId() > 0)
        {
            new_location[0] = cell_location[0];
        }

        // Move node to new location
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);
        p_node->rGetModifiableLocation() = new_location;
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
