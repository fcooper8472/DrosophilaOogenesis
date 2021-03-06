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

#ifndef GERMARIUMBOUNDARYCONDITION_HPP_
#define GERMARIUMBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A cell population boundary condition class, which restricts nodes to lie
 * on the surface of a set of crypts and a villus.
 */
class GermariumBoundaryCondition : public AbstractCellPopulationBoundaryCondition<3>
{
private:

    /** The radius of the crypt bases. */
    double mRadiusOfGermarium;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<3> >(*this);
        archive & mRadiusOfGermarium;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param germariumRadius the radius of the germarium
     */
    GermariumBoundaryCondition(AbstractCellPopulation<3>* pCellPopulation, double germariumRadius);

    /**
     * @return mRadiusOfGermarium the radius of the Germarium
     */
    double GetRadiusOfGermarium() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(GermariumBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a GermariumBoundaryCondition.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const GermariumBoundaryCondition* t, const unsigned int file_version)
{
    // Get data required to construct instance
    const AbstractCellPopulation<3>* const p_cell_population = t->GetCellPopulation();
    double germarium_radius = t->GetRadiusOfGermarium();

    // Archive
    ar << p_cell_population;
    ar << germarium_radius;
}

/**
 * De-serialize constructor parameters and initialize a GermariumBoundaryCondition.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, GermariumBoundaryCondition* t, const unsigned int file_version)
{
    // Make variables required to construct new instance with constructor
    AbstractCellPopulation<3>* p_cell_population;
    double germarium_radius;

    // Populate these variables
    ar >> p_cell_population;
    ar >> germarium_radius;

    // Invoke constructor to initialise instance
    ::new(t)GermariumBoundaryCondition(p_cell_population, germarium_radius);
}
}
} // namespace ...

#endif /*GERMARIUMBOUNDARYCONDITION_HPP_*/
