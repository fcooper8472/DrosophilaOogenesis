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

#ifndef DIFFERENTIALADHESIONGERMARIUMFORCE_HPP_
#define DIFFERENTIALADHESIONGERMARIUMFORCE_HPP_

#include "GeneralisedLinearSpringForce.hpp"

/**
 * A class for a simple two-body differential adhesion force law between
 * labelled and unlabelled cells (as defined by the CellLabel cell
 * property).
 *
 * Designed for use in node and mesh-based simulations.
 *
 * \todo #2266 - throw exceptions if using other cell population objects?
 * \todo #2266 - override CalculateForceBetweenNodes() to use a default rest length of 1.0 for all springs?
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class DifferentialAdhesionGermariumForce : public GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>
{
private :

    /**
     * Multiplier for interactions between aggregate cells.
     */
    double mAggregateAggregateSpringConstantMultiplier;

    /**
     * Multiplier for interactions between Follicle cells.
     */
    double mFollicleFollicleSpringConstantMultiplier;

    /**
     * Multiplier for interactions between Follicle and aggregate cells.
     */
    double mFollicleAggregateSpringConstantMultiplier;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mAggregateAggregateSpringConstantMultiplier;
        archive & mFollicleFollicleSpringConstantMultiplier;
        archive & mFollicleAggregateSpringConstantMultiplier;
    }

public :

    /**
     * Constructor.
     */
    DifferentialAdhesionGermariumForce();

    /**
     * Overridden VariableSpringConstantMultiplicationFactor() method.
     *
     * This method takes account of the distinct spring constants present
     * for homotypic (labelled-labelled and unlabelled-unlabelled) and
     * heterotypic (labelled-unlabelled) interactions between neighbouring
     * cells.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                      unsigned nodeBGlobalIndex,
                                                      AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                      bool isCloserThanRestLength);

    /**
     * All the set and get methods for the different spring constant multipliers in this class
     */
    double GetAggregateAggregateSpringConstantMultiplier();

    double GetFollicleFollicleSpringConstantMultiplier();

    double GetFollicleAggregateSpringConstantMultiplier();

    void SetAggregateAggregateSpringConstantMultiplier(double aggregateAggregateSpringConstantMultiplier);

    void SetFollicleFollicleSpringConstantMultiplier(double follicleFollicleSpringConstantMultiplier);

    void SetFollicleAggregateSpringConstantMultiplier(double follicleAggregateSpringConstantMultiplier);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DifferentialAdhesionGermariumForce)

#endif /*DIFFERENTIALADHESIONGERMARIUMFORCE_HPP_*/
