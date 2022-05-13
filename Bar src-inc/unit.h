#ifndef UNIT_H_
#define UNIT_H_

//
// Define unit constants
//											Basic units (compatible with USPC7100)
#define UNIT_NONE							-1
#define UNIT_us								0
#define UNIT_mm								1
#define UNIT_in								2
//											Time units
#define UNIT_ns								100
#define UNIT_ms								101
#define UNIT_s								102
//											Distance units
#define UNIT_um								200
#define UNIT_m								201
//											Frequency units
#define UNIT_Hz								300 
#define UNIT_kHz							301
#define UNIT_MHz							302
//											Velocity units
#define UNIT_mps							400
#define UNIT_inpus							401
//											Angle units
#define UNIT_deg							500
#define UNIT_rad							501
//											Temperture units
#define UNIT_Celsius						600
#define UNIT_Fahrenheit						601

class Unit
{
public:
	static bool ChangeUnit(double* Value, int FromUnit, int ToUnit);
	static double ChangeUnit(double Value, int FromUnit, int ToUnit);
	static int GetUnitFromString(char* UnitString);
	static double Convert(double value, int FromUnit, int ToUnit, double velocity);
};

#endif /* UNIT_H_ */