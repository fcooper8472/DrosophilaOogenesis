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

#ifndef MULTIPLECRYPTGEOMETRYBOUNDARYCONDITION_HPP_
#define MULTIPLECRYPTGEOMETRYBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A cell population boundary condition class, which restricts nodes to lie
 * on the surface of a set of crypts and a villus.
 */
class MultipleCryptGeometryBoundaryCondition : public AbstractCellPopulationBoundaryCondition<3>
{
private:

    /** The radius of the crypt bases. */
    double mRadiusOfCrypt;

    /** The length of the crypts. */
    double mLengthOfCrypt;

    /** The radius of the villus base. */
    double mRadiusOfVillus;

    /** The length of the villus. */
    double mLengthOfVillus;

    /** The domain width in x and y directions */
    double mDomainWidth;

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
        archive & mRadiusOfCrypt;
        archive & mLengthOfCrypt;
        archive & mRadiusOfVillus;
        archive & mLengthOfVillus;
        archive & mDomainWidth;
    }

    /**
     * Move a cell onto a crypt
     *
     * @param rCellLocation  the old cell location after a movement
     * @param rBaseCentre  the location of the crypt base
     * @return the new cell location after BC enforced
     */
    c_vector<double, 3> MoveToCrypt(const c_vector<double,3>& rCellLocation, const c_vector<double,3>& rBaseCentre);


public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param cryptRadius the radius of the crypt bases (and top rim)
     * @param cryptLength the length of the crypts
     * @param villusRadius the radius of the villus top (and base rim)
     * @param villusLength the length of the villus
     * @param domainWidth the width (in x and y directions) of the domain.
     */
    MultipleCryptGeometryBoundaryCondition(AbstractCellPopulation<3>* pCellPopulation,
                                    double cryptRadius,
                                    double cryptLength,
                                    double villusRadius,
                                    double villusLength,
                                    double domainWidth);

    /**
     * @return The radii of the crypts.
     */
    double GetRadiusOfCrypt() const;

    /**
     * @return The length of the crypts.
     */
    double GetLengthOfCrypt() const;

    /**
     * @return The radius of the villus.
     */
    double GetRadiusOfVillus() const;

    /**
     * @return The length of the villus.
     */
    double GetLengthOfVillus() const;

    /**
     * @return The width of the domain (in x and y directions).
     */
    double GetDomainWidth() const;

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
CHASTE_CLASS_EXPORT(MultipleCryptGeometryBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MultipleCryptGeometryBoundaryCondition.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const MultipleCryptGeometryBoundaryCondition* t, const unsigned int file_version)
{
    // Get data required to construct instance
    const AbstractCellPopulation<3>* const p_cell_population = t->GetCellPopulation();
    double crypt_radius = t->GetRadiusOfCrypt();
    double crypt_length = t->GetLengthOfCrypt();
    double villus_radius = t->GetRadiusOfVillus();
    double villus_length = t->GetLengthOfVillus();
    double domain_width = t->GetDomainWidth();

    // Archive
    ar << p_cell_population;
    ar << crypt_radius;
    ar << crypt_length;
    ar << villus_radius;
    ar << villus_length;
    ar << domain_width;
}

/**
 * De-serialize constructor parameters and initialize a MultipleCryptGeometryBoundaryCondition.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, MultipleCryptGeometryBoundaryCondition* t, const unsigned int file_version)
{
    // Make variables required to construct new instance with constructor
    AbstractCellPopulation<3>* p_cell_population;
    double crypt_radius;
    double crypt_length;
    double villus_radius;
    double villus_length;
    double domain_width;

    // Populate these variables
    ar >> p_cell_population;
    ar >> crypt_radius;
    ar >> crypt_length;
    ar >> villus_radius;
    ar >> villus_length;
    ar >> domain_width;

    // Invoke constructor to initialise instance
    ::new(t)MultipleCryptGeometryBoundaryCondition(p_cell_population, crypt_radius, crypt_length,
                                                   villus_radius, villus_length, domain_width);
}
}
} // namespace ...

#endif /*MULTIPLECRYPTGEOMETRYBOUNDARYCONDITION_HPP_*/
