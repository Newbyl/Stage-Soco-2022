#ifndef CBAR_H__
#define CBAR_H__

#define _CRT_SECURE_NO_WARNINGS

#define M_PI 3.14159265358979323846264338327950288
#define _OPTIMIZATION

#include "IPlugin.h"
#include "framework.h"
#include "unit.h"
#include "error_codes.h"

#include <math.h>
#include <vector>


class Cbar : public IPlugin
{
public:
    Cbar();
    ~Cbar();

    int Open();
    int Close();
    int Set(const char* param_name, int unit, int* value);
    int Set(const char* param_name, int unit, double* value);
    int Set(const char* param_name, int unit, char* value);
    int Set(const char* param_name, int unit, int* dim1, int* value);
    int Set(const char* param_name, int unit, int* dim1, double* value);
    int Set(const char* param_name, int unit, int* dim1, int* dim2, int* value);
    int Set(const char* param_name, int unit, int* dim1, int* dim2, double* value);
    int Set(const char* param_name, int unit, int* dim1, int* dim2, int* dim3, int* value);
    int Set(const char* param_name, int unit, int* dim1, int* dim2, int* dim3, double* value);
    int Get(const char* param_name, int unit, int* value);
    int Get(const char* param_name, int unit, double* value);
    int Get(const char* param_name, int unit, int* dim1, int* value);
    int Get(const char* param_name, int unit, int* dim1, double* value);
    int Get(const char* param_name, int unit, int* dim1, char* value);
    int Get(const char* param_name, int unit, int* dim1, int* dim2, int* value);
    int Get(const char* param_name, int unit, int* dim1, int* dim2, double* value);
    int Get(const char* param_name, int unit, int* dim1, int* dim2, int* dim3, int* value);
    int Get(const char* param_name, int unit, int* dim1, int* dim2, int* dim3, double* value);
    int ExecSync(const char* action);
    int ExecAsync(const char* action);

private:
    bool opened;                    // Indicates if this class is ready to use or not.
    bool calculationDone;           // Indicates if the calculation is already done.
    int numberOfElements;           // Number of elements to calculate the laws.
    
    double resolution;              // The resolution for interface points, default = 0.1mm.
    double diameter;                // Diameter of the bar in millimeters.

    double centerApertureX;         // x,y,z coordinates of the center of the aperture.
    double centerApertureY;         //
    double centerApertureZ;         //

    struct ELEMENTS {               // Spacial coordinates of the elements.
        struct COORDINATES {
            double* x;              // in millimeters
            double* y;              //
            double* z;              //
        }coordinates;
    }elements;

    struct MATERIAL {               // Material characteristic 
        double velocity;            // in m/s
    }material;

    struct COUPLING {               // Coupling (interface) characteristic
        double velocity;            // in m/s
        double height;              // in millimeters
    }coupling;

    enum class ANGLE_TYPE {         // Angle type of notches (0 = INCIDENT, 1 = TRANSMITED)
        INCIDENT,
        TRANSMITED
    }angleType;

    enum class DEFECT_TYPE {         // Type of defect you want to see (0 = NOTCHE, 1 = FBH)
        NOTCHE,
        FBH
    }defectType;

    int numberOfTargets;            // Number of targets determines the number of laws of positions of remarkable points to calculate  

    typedef struct _TARGETS {
        double* tilts;              // Defines the tilt angle for each target. In degree.
        double* notchesAngles;      // Defines the notche angle for each target. In degree.
        double* positions;          // Defines the relative position for each target. In millimeters. Depends on positionType paremeter.
    }TARGET, *PTARGET;
    TARGET targets;                 // List of targets whose laws and remarkable points of their routes must be calculated.

    typedef struct _LAW {
        double* delays;             // An array of delays calculated by this class to apply on elements to obtain an acoustic answer on a target. In usec.
    }LAW, * PLAW;
    PLAW laws;                      // An array of laws calculated by this class to apply to obtain an acoustic answer on targets.

    typedef struct _PATH {
        double x[3];                // Array of X coordinates of remarkable points (source, interface, default) for a given target. In millimeters.
        double y[3];                // Array of Y coordinates of remarkable points (source, interface, default) for a given target. In millimeters.
        double z[3];                // Array of Z coordinates of remarkable points (source, interface, default) for a given target. In millimeters.
    }PATH, * PPATH;                 
    PPATH paths;                    // An array of remarkable points (source, interface, defect) for each target

    int Calculate();                // Calculates laws and paths
    double maxArray(double *array, int size);   // return the max element of an array
    double minArray(double *array, int size);   // return the max element of an array
    std::vector<double*> fbhBuilder(double barDiameter2);
    std::vector<double*> notcheBuilder(double barDiameter2);
};


#if defined(_WINDOWS)
    //  Microsoft
    #define EXPORT extern "C" __declspec(dllexport)
    #define IMPORT extern "C" __declspec(dllimport)
#elif defined(_LINUX)
    //  GCC
    #define EXPORT extern "C" __attribute__((visibility("default")))
    #define IMPORT
#else
    //  do nothing and hope for the best?
    #define EXPORT
    #define IMPORT
    #pragma warning Unknown dynamic link import/export semantics.
#endif

// Finally export a creation and deletion function (which you can get using LoadLibrary/GetProcAddress)
EXPORT IPlugin* CreatePluginInstance()
{
    return new Cbar();
}
EXPORT void ReleasePluginInstance(IPlugin* p)
{
    p->Close();
}

#endif /* CBAR_H__ */
