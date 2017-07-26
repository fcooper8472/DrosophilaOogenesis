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

#include "GermariumDivisionRule.hpp"

#include "CellLabel.hpp"
#include "DrosophilaOogenesisEnumerations.hpp"

#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GermariumDivisionRule<ELEMENT_DIM, SPACE_DIM>::GermariumDivisionRule(c_vector<double, SPACE_DIM>& rAggregateOffset,
                                                                     c_vector<double, SPACE_DIM>& rFollicleOffset)
        : mAggregateOffset(rAggregateOffset),
          mFollicleOffset(rFollicleOffset)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& GermariumDivisionRule<ELEMENT_DIM,SPACE_DIM>::rGetAggregateOffset() const
{
    return mAggregateOffset;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& GermariumDivisionRule<ELEMENT_DIM,SPACE_DIM>::rGetFollicleOffset() const
{
    return mFollicleOffset;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > GermariumDivisionRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    c_vector<double, SPACE_DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell);

    // We get the cell label from the current cell and obtain its colour
    auto abstract_property = pParentCell->rGetCellPropertyCollection().GetProperties<CellLabel>().GetProperty();
    unsigned cell_colour = boost::static_pointer_cast<CellLabel>(abstract_property)->GetColour();

    c_vector<double, SPACE_DIM> daughter_position;

    if (cell_colour == TYPE_AGGREGATE)
    {
        daughter_position = parent_position + mAggregateOffset;
        MARK;
        PRINT_VECTOR(daughter_position);
    }
    else if (cell_colour == TYPE_FOLLICLE)
    {
        daughter_position = parent_position + mFollicleOffset;
        MARK;
        PRINT_VECTOR(daughter_position);
    }
    else
    {
        NEVER_REACHED;
    }

    return std::make_pair(parent_position, daughter_position);
}

// Explicit instantiation
template class GermariumDivisionRule<1,1>;
template class GermariumDivisionRule<1,2>;
template class GermariumDivisionRule<2,2>;
template class GermariumDivisionRule<1,3>;
template class GermariumDivisionRule<2,3>;
template class GermariumDivisionRule<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(GermariumDivisionRule)