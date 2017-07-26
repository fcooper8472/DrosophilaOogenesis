/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "DifferentialAdhesionGermariumForce.hpp"
#include "DrosophilaOogenesisEnumerations.hpp"
#include "CellLabel.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::DifferentialAdhesionGermariumForce()
   : GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(),
     mAggregateAggregateSpringConstantMultiplier(0.0),
     mFollicleFollicleSpringConstantMultiplier(1.0),
     mFollicleAggregateSpringConstantMultiplier(2.0)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
    bool isCloserThanRestLength)
{

    if (isCloserThanRestLength)
    {
        return 1.0;
    }
    else
    {
        // Determine which (if any) of the cells corresponding to these nodes are labelled

        CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
        auto abstract_property_a = p_cell_A->
                rGetCellPropertyCollection().GetProperties<CellLabel>().GetProperty();
        unsigned cell_colour_a = boost::static_pointer_cast<CellLabel>(abstract_property_a)->GetColour();


        CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);
        auto abstract_property_b = p_cell_B->
                rGetCellPropertyCollection().GetProperties<CellLabel>().GetProperty();
        unsigned cell_colour_b = boost::static_pointer_cast<CellLabel>(abstract_property_b)->GetColour();

        // For heterotypic interactions, scale the spring constant by mHeterotypicSpringConstantMultiplier
        if (cell_colour_a == TYPE_AGGREGATE && cell_colour_b == TYPE_AGGREGATE)
        {
            return mAggregateAggregateSpringConstantMultiplier;
        }
        else if (cell_colour_a == TYPE_FOLLICLE && cell_colour_b == TYPE_FOLLICLE)
        {
            return mFollicleFollicleSpringConstantMultiplier;
        }
        else if ( cell_colour_a != cell_colour_b )
        {
            return mFollicleAggregateSpringConstantMultiplier;
        }
        else
        {
            NEVER_REACHED;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::GetAggregateAggregateSpringConstantMultiplier()
{
    return mAggregateAggregateSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::GetFollicleFollicleSpringConstantMultiplier()
{
    return mFollicleFollicleSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::GetFollicleAggregateSpringConstantMultiplier()
{
    return mFollicleAggregateSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::SetAggregateAggregateSpringConstantMultiplier(
        double aggregateAggregateSpringConstantMultiplier)
{
    assert(aggregateAggregateSpringConstantMultiplier > 0.0);
    mAggregateAggregateSpringConstantMultiplier = aggregateAggregateSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::SetFollicleFollicleSpringConstantMultiplier(
        double follicleFollicleSpringConstantMultiplier)
{
    assert(follicleFollicleSpringConstantMultiplier > 0.0);
    mFollicleFollicleSpringConstantMultiplier = follicleFollicleSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::SetFollicleAggregateSpringConstantMultiplier(
        double follicleAggregateSpringConstantMultiplier)
{
    assert(follicleAggregateSpringConstantMultiplier > 0.0);
    mFollicleAggregateSpringConstantMultiplier = follicleAggregateSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialAdhesionGermariumForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AggregateAggregateSpringConstantMultiplier>" << mAggregateAggregateSpringConstantMultiplier << "</AggregateAggregateSpringConstantMultiplier>\n";
    *rParamsFile << "\t\t\t<FollicleFollicleSpringConstantMultiplier>" << mFollicleFollicleSpringConstantMultiplier << "</FollicleFollicleSpringConstantMultiplier>\n";
    *rParamsFile << "\t\t\t<FollicleAggregateSpringConstantMultiplier>" << mFollicleAggregateSpringConstantMultiplier << "</FollicleAggregateSpringConstantMultiplier>\n";

    // Call direct parent class
    GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class DifferentialAdhesionGermariumForce<1,1>;
template class DifferentialAdhesionGermariumForce<1,2>;
template class DifferentialAdhesionGermariumForce<2,2>;
template class DifferentialAdhesionGermariumForce<1,3>;
template class DifferentialAdhesionGermariumForce<2,3>;
template class DifferentialAdhesionGermariumForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DifferentialAdhesionGermariumForce)
