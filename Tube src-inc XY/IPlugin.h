#ifndef IPLUGIN_H_
#define IPLUGIN_H_

#define PLUGIN_NO_ERROR						 0
#define PLUGIN_EXCEPTION					-1

#define PLUGIN_NULL_POINTER					-10
#define PLUGIN_INVALID_HANDLE				-11
#define PLUGIN_MEMORY_ALLOCATION_FAILED		-12

#define PLUGIN_NOT_OPEN						-100
#define PLUGIN_ALREADY_OPEN					-101
#define PLUGIN_NOT_FOUND                    -102
#define PLUGIN_LICENSE_EXPIRED              -103

#define	PLUGIN_INVALID_ARGUMENT				-110
#define PLUGIN_UNKNOWN_PARAMETER			-111

class IPlugin
{
public:
   //virtual ~IPlugin() = 0;

    virtual int Open() = 0;
    virtual int Close() = 0;
    
    virtual int Set(const char* param_name, int unit, int* value) = 0;
    virtual int Set(const char* param_name, int unit, double* value) = 0;
    virtual int Set(const char* param_name, int unit, char* value) = 0;
    virtual int Set(const char* param_name, int unit, int *dim1, int* value) = 0;
    virtual int Set(const char* param_name, int unit, int *dim1, double* value) = 0;
    virtual int Set(const char* param_name, int unit, int *dim1, int *dim2, int* value) = 0;
    virtual int Set(const char* param_name, int unit, int *dim1, int *dim2, double* value) = 0;
    virtual int Set(const char* param_name, int unit, int *dim1, int *dim2, int *dim3, int* value) = 0;
    virtual int Set(const char* param_name, int unit, int *dim1, int *dim2, int *dim3, double* value) = 0;

    virtual int Get(const char* param_name, int unit, int* value) = 0;
    virtual int Get(const char* param_name, int unit, double* value) = 0;
    virtual int Get(const char* param_name, int unit, int *dim1, int* value) = 0;
    virtual int Get(const char* param_name, int unit, int* dim1, double* value) = 0;
    virtual int Get(const char* param_name, int unit, int* dim1, char* value) = 0;
    virtual int Get(const char* param_name, int unit, int *dim1, int *dim2, int* value) = 0;
    virtual int Get(const char* param_name, int unit, int *dim1, int *dim2, double* value) = 0;
    virtual int Get(const char* param_name, int unit, int *dim1, int *dim2, int *dim3, int* value) = 0;
    virtual int Get(const char* param_name, int unit, int *dim1, int *dim2, int *dim3, double* value) = 0;

    virtual int ExecSync(const char* action) = 0;
    virtual int ExecAsync(const char* action) = 0;
    ~IPlugin();
};

#endif /* IPLUGIN_H_ */

