#ifndef __SL2CFOAM_DBG_H__
#define __SL2CFOAM_DBG_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

/**********************************************************************/

////////////////////////////////////////////////////////////////////////
// Functions useful for debugging.
////////////////////////////////////////////////////////////////////////

#define VARNAME(x) #x

// Print the name and values of variables.
#define DPI(x) printf("DEBUG - "VARNAME(x)": %i\n", (int)x);
#define DPD(x) printf("DEBUG - "VARNAME(x)": %g\n", (double)x);
#define DPS(x) printf("DEBUG - "VARNAME(x)": %s\n", (char*)x);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_DBG_H__*/