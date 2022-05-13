#ifndef CPLATE_H__
#define CPLATE_H__

#include "IPlugin.h"

class Cplate : public IPlugin
{
public:
    Cplate();
    ~Cplate();

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
    
    int probeType;
    double resolution;              // The resolution for interface points, default = 0.1mm
    double pitch;
    

    struct ELEMENTS {               // Spacial coordinates of the elements.
        struct COORDINATES {
            double* x;              // in millimeters
            double* y;              //
            double* z;              //
            double* i;
        }coordinates;
    }elements;

    struct MATERIAL {               // Material characteristic 
        double velocity;            // in m/s
    }material;

    struct COUPLING {               // Coupling (interface) characteristic
        double velocity;            // in m/s
        double height;              // in millimeters
    }coupling;

    enum class FLOW_POSITION_TYPE { // How is defined the position of flow
        PATH,                       // Defined on ultrasonic path
        DEPTH                       // Defined from the surface of the material
    }positionType;

    int numberOfTargets;            // Number of targets determines the number of laws of positions of remarkable points to calculate  

    typedef struct _TARGETS {
        double* tilts;              // Defines the tilt angle for each target. In degree.
        double* skews;              // Defines the skew angle for each target. In degree.
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
    double maxArray(const double *array, int size);
    double minArray(const double *array, int size);
    double *append(double *ar1, double *ar2, int len1, int len2);
    int minArrayIndex(const double *array, int size);
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


#endif /* CPLATE_H__ */
