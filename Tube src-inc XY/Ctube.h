#ifndef CTUBE_H__
#define CTUBE_H__

#include "IPlugin.h"
#include "unit.h"
#include "error_codes.h"
#include "framework.h"

#include <vector>
#include <cmath>
#include <iostream>

// We include Eigen headers here, which is a library for linear algebra,
// we use it here to manipulate matrix easily.
// Eigen folder with all header files is in the VS project in "inc" folder.
#include "Eigen/LU"
#include "Eigen/Core"

// Definition of PI used in the .cpp file
#define M_PI 3.14159265358979323846

# define _OPTIMIZATION

class Ctube : public IPlugin
{
public:
    Ctube();
    ~Ctube();

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
    double resolution;

    double diameter;                // Diameter of the tube.
    double thickness;               // Thickness of the tube.

    double centerApertureX;
    double centerApertureY;
    double centerApertureZ;

    enum class DEFECT_TYPE {         // Type of defect you want to see (0 = NOTCHE, 1 = FBH)
        NOTCHE,
        FBH
    }defectType;

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

    struct FOCAL {                  // Complete path (also in water)
        struct LENGTH {
            double coupling;        // in millimeters
        }length;
    }focal;


    int numberOfTargets;            // Number of targets determines the number of laws of positions of remarkable points to calculate  

    typedef struct _TARGETS {
        double* tilts;              // Defines the tilt angle for each target. In degree.
        double* skews;              // Defines the skew angle for each target. In degree.
        double* positions;
    }TARGET, * PTARGET;
    TARGET targets;                 // List of targets whose laws and remarkable points of their routes must be calculated.

    typedef struct _LAW {
        double* delays;             // An array of delays calculated by this class to apply on elements to obtain an acoustic answer on a target. In usec.
    }LAW, * PLAW;
    PLAW laws;                      // An array of laws calculated by this class to apply to obtain an acoustic answer on targets.

    typedef struct _PATH {
        double x[4];                // Array of X coordinates of remarkable points (source, interface, default) for a given target. In millimeters.
        double y[4];                // Array of Y coordinates of remarkable points (source, interface, default) for a given target. In millimeters.
        double z[4];                // Array of Z coordinates of remarkable points (source, interface, default) for a given target. In millimeters.
    }PATH, * PPATH;
    PPATH paths;                    // An array of remarkable points (source, interface, defect) for each target


    double* alphaS;                 // Deflection angle.
    double* xi3D;                   // Coordinates of the interface on X axis.
    double* yi3D;                   // Coordinates of the interface on Y axis.
    double* zi3D;                   // Coordinates of the interface on Z axis.

    int Calculate();                                                    // Calculates laws and paths
    double maxArray(const double* array, int size);                     // Find the max element of an array
    double minArray(const double* array, int size);                     // Find the min element of an array
    std::vector<double> newElipse(double skew, double alphaI);          // à demander
    std::vector<double*> fbhBuilder(double barDiameter2);
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

//Finally export a creation and deletion function (which you can get using LoadLibrary/GetProcAddress)
EXPORT IPlugin* CreatePluginInstance()
{
    return new Ctube();
}
EXPORT void ReleasePluginInstance(IPlugin* p)
{
    p->Close();
}

#endif /* CTUBE_H__ */
