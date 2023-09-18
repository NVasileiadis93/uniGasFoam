#include "uspCloud.H"
#include "CloudFunctionObject.H"
//#include "FacePostProcessing.H"
#include "uspFaceTrackerFO.H"

#include "runTimeSelectionTables.H"
/*
    defineNamedTemplateTypeNameAndDebug
    (
        Foam::CloudFunctionObject<Foam::uspCloud>,
        0
    );
    namespace Foam
    {
        defineTemplateRunTimeSelectionTable
        (
            CloudFunctionObject<uspCloud>,
            dictionary
        );
    }
*/
makeCloudFunctionObject(uspCloud);

namespace Foam
{
//makeCloudFunctionObjectType(FacePostProcessing, uspCloud);
//makeCloudFunctionObjectTypeNT(uspFaceTrackerFO, uspCloud);
}

