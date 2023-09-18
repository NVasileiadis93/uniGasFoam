#ifndef CloudWithModels_H
#define CloudWithModels_H

#include "Cloud.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "cloudSolution.H"

namespace Foam
{
template<class Type>
class CloudWithModels
:
    public Cloud<Type>
{
public:

    // Public typedefs

        //- Type of particle the cloud was instantiated for
        typedef Type parcelType;

        typedef CloudWithModels<Type> baseCloudType;


protected:

    // Protected data

        //- Reference to the fvMesh
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Dictionary of output properties
        IOdictionary outputProperties_;

        //- Solution properties
        cloudSolution solution_;

        //- Cell occupancy information for each parcel, (demand driven)
        autoPtr<List<DynamicList<Type*>>> cellOccupancyPtr_;


    // Protected Member Functions


public:

    // Constructors

        //- Construct from mesh and a list of particles
        CloudWithModels
        (
            const fvMesh& mesh,
            const word& cloudName,
            const IDLList<Type>& particles
        );

        //- Construct from mesh by reading from file with given cloud instance
        //  Optionally disable checking of class name for post-processing
        CloudWithModels
        (
            const fvMesh& pMesh,
            const word& cloudName,
            const bool checkClass = true
        );


    // Member Functions

        //- Return a reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return particle properties dictionary
        inline const IOdictionary& particleProperties() const;

        //- Return output properties dictionary
        inline const IOdictionary& outputProperties() const;

        //- Return non-const access to the output properties dictionary
        inline IOdictionary& outputProperties();

        //- Return const access to the solution properties
        inline const cloudSolution& solution() const;


    // Cell occupancy

        //- Return the cell occupancy information for each parcel, non-const
        //- access, the caller is responsible for updating it for its own
        //- purposes if particles are removed or created.
        inline List<DynamicList<parcelType*>>& cellOccupancy();

        //- Build the cellOccupancy
        void buildCellOccupancy();

        //- Update/build the cellOccupancy if it has already been used
        void updateCellOccupancy();

        //- Insert parcel
        void insertParcelInCellOccupancy(parcelType* p);

        //- Remove parcel
        void removeParcelFromCellOccupancy
        (
            const label parceli,
            const label celli
        );
};

} // End namespace Foam

#ifdef NoRepository
    #include "CloudWithModels.C"
#endif

#include "CloudWithModelsI.H"

#endif
