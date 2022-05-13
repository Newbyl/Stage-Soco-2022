#include "framework.h"
#include "unit.h"
#include <math.h>
#include <string.h>
#include "framework.h"
#define __min(a,b) (((a) < (b)) ? (a) : (b))

#define MY_PI 3.14159265358979323846

bool Unit::ChangeUnit(double* Value, int FromUnit, int ToUnit)
{
	float from, to;

	// Is it necessary to convert?
	if (FromUnit == ToUnit)
		// No, it is not necessary
		return true;

	switch (FromUnit)
	{
	case UNIT_ns:
		from = 1e-9;
		break;
	case UNIT_us:
		//from = 1e-6;
		from = 0.000001;
		break;
	case UNIT_ms:
		from = 1e-3;
		break;
	case UNIT_s:
		from = 1.0;
		break;
	case UNIT_um:
		from = 1e-6;
		break;
	case UNIT_mm:
		from = 1e-3;
		break;
	case UNIT_m:
		from = 1.0;
		break;
	case UNIT_in:
		from = 25.4e-3;
		break;
	case UNIT_Hz:
		from = 1.0;
		break;
	case UNIT_kHz:
		from = 1e3;
		break;
	case UNIT_MHz:
		from = 1e6;
		break;
	case UNIT_mps:
		from = 1.0;
		break;
	case UNIT_inpus:
		from = 25400.0; // m/s
		break;
	case UNIT_deg:
		from = 1.0;
		break;
	case UNIT_rad:
		from = 180.0 / MY_PI;
		break;
	case UNIT_Celsius:
		from = 33.8;
		break;
	case UNIT_Fahrenheit:
		from = 1.0;
		break;
	default:
		return false;
	}
	switch (ToUnit)
	{
	case UNIT_ns:
		to = 1e-9;
		break;
	case UNIT_us:
		to = 1e-6;
		break;
	case UNIT_ms:
		to = 1e-3;
		break;
	case UNIT_s:
		to = 1.0;
		break;
	case UNIT_um:
		to = 1e-6;
		break;
	case UNIT_mm:
		to = 1e-3;
		break;
	case UNIT_m:
		to = 1.0;
		break;
	case UNIT_in:
		to = 25.4e-3;
		break;
	case UNIT_Hz:
		to = 1.0;
		break;
	case UNIT_kHz:
		to = 1e3;
		break;
	case UNIT_MHz:
		to = 1e6;
		break;
	case UNIT_mps:
		to = 1.0;
		break;
	case UNIT_inpus:
		to = 25400.0;
		break;
	case UNIT_deg:
		to = 1.0;
		break;
	case UNIT_rad:
		to = 180.0 / MY_PI;
		break;
	case UNIT_Celsius:
		if (FromUnit == UNIT_Fahrenheit)
		{
			*Value = ((*Value) - 32.0) * 5.0 / 9.0;
			*Value = (double)((float)(*Value));
			return true;
		}
		return false;
		break;
	case UNIT_Fahrenheit:
		if (FromUnit == UNIT_Celsius)
		{
			*Value = 9.0 / 5.0 * (*Value) + 32.0;
			*Value = (double)((float)(*Value));
			return true;
		}
		return false;
		break;
	default:
		return false;
	}
	float conversion_factor = from / to;
	*Value = *Value * conversion_factor;
	//*Value = (double)((float)(*Value));
	return true;
}

double Unit::ChangeUnit(double Value, int FromUnit, int ToUnit)
{
	ChangeUnit(&Value, FromUnit, ToUnit);
	return Value;
}

int Unit::GetUnitFromString(char* UnitString)
{
	char Unit[20];

	if (UnitString == NULL) return UNIT_NONE;

	// Delete first spaces
	int i = 0;
	while (UnitString[i] == ' ' && i < 20)i++;
	if (i < 20)
	{
		memset(Unit, 20, 0);
		strncpy(Unit, UnitString + i, __min(strlen(UnitString + i), 19));
	}
	else
		return UNIT_NONE;

	// Comparaison is ordered by most popular units to optimize the research time

	if (strcmp(Unit, "us") == 0) return UNIT_us;
	if (strcmp(Unit, "mm") == 0) return UNIT_mm;
	if (strcmp(Unit, "in") == 0) return UNIT_in;
	if (strcmp(Unit, "m/s") == 0) return UNIT_mps;
	if (strcmp(Unit, "in/us") == 0) return UNIT_inpus;

	if (strcmp(Unit, "ns") == 0) return UNIT_ns;
	if (strcmp(Unit, "ms") == 0) return UNIT_ms;
	if (strcmp(Unit, "s") == 0) return UNIT_s;

	if (strcmp(Unit, "um") == 0) return UNIT_um;
	if (strcmp(Unit, "m") == 0) return UNIT_m;

	if (strcmp(Unit, "Hz") == 0) return UNIT_Hz;
	if (strcmp(Unit, "kHz") == 0) return UNIT_kHz;
	if (strcmp(Unit, "MHz") == 0) return UNIT_MHz;

	if (strcmp(Unit, "s^-1 m") == 0) return UNIT_mps;
	if (strcmp(Unit, "us^-1 in") == 0) return UNIT_inpus;

	if (strcmp(Unit, "deg") == 0) return UNIT_deg;
	if (strcmp(Unit, "°") == 0) return UNIT_deg;
	if (strcmp(Unit, "rad") == 0) return UNIT_rad;

	return UNIT_NONE;
}

double Unit::Convert(double value, int FromUnit, int ToUnit, double velocity)										// Note: velocity must be m/s
{
	bool FromTimeToDistance;

	if ((FromUnit == UNIT_us || FromUnit == UNIT_ns || FromUnit == UNIT_ms || FromUnit == UNIT_s) &&
		(ToUnit == UNIT_mm || ToUnit == UNIT_in))
	{	// Convert time to distance
		ChangeUnit(&value, FromUnit, UNIT_us);																		// Change time to us
		switch (ToUnit)
		{
		case UNIT_mm:
			value = (value * (velocity / 2.0)) / 1000.0;															// Convert us to mm
			break;
		case UNIT_in:
			value = ((value * (velocity / 2.0)) / 1000.0) / 25.4;													// Convert us to inch
			break;
		}
		return value;
	}
	if ((ToUnit == UNIT_us || ToUnit == UNIT_ns || ToUnit == UNIT_ms || ToUnit == UNIT_s) &&
		(FromUnit == UNIT_mm || FromUnit == UNIT_in))
	{	// Convert distance to time
		ChangeUnit(&value, FromUnit, UNIT_mm);																		// Change distance to mm
		value = (value * 1000.0) / (velocity / 2.0);																// Convert mm to us
		switch (ToUnit)
		{
		case UNIT_ns:
			ChangeUnit(&value, UNIT_us, UNIT_ns);																	// Change us to ns
			break;
		case UNIT_ms:
			ChangeUnit(&value, UNIT_us, UNIT_ms);																	// Change us to ms
			break;
		case UNIT_s:
			ChangeUnit(&value, UNIT_us, UNIT_s);																	// Change us to s
			break;
		}
		return value;
	}
	if ((FromUnit == UNIT_us || FromUnit == UNIT_ns || FromUnit == UNIT_ms || FromUnit == UNIT_s) &&
		(ToUnit == UNIT_us || ToUnit == UNIT_ns || ToUnit == UNIT_ms || ToUnit == UNIT_s))
	{	// Convert time to time
		ChangeUnit(&value, FromUnit, ToUnit);																		// Change time unit
		return value;
	}
	if ((FromUnit == UNIT_mm || FromUnit == UNIT_um || FromUnit == UNIT_m || FromUnit == UNIT_in) &&
		(ToUnit == UNIT_mm || ToUnit == UNIT_um || ToUnit == UNIT_m || ToUnit == UNIT_in))
	{	// Convert distance to distance
		ChangeUnit(&value, FromUnit, ToUnit);																		// Change time unit
		return value;
	}

	// No conversion
	return value;
}
