#ifndef __SL2CFOAM_ERROR_H__
#define __SL2CFOAM_ERROR_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////
// Functions to handle exceptions.
////////////////////////////////////////////////////////////////////////

// Prints the file and line of the error location, an optional message
// and then interrupts program.
#define error(format, ...) {\
        fprintf(stderr, "sl2cfoam: ERROR at file %s, line %d. Terminating...\n", __FILE__, __LINE__);\
        fprintf(stderr, "        : ");\
        fprintf(stderr, format, ##__VA_ARGS__);\
        fprintf(stderr, "\n");\
        exit(EXIT_FAILURE);\
        }

// Prints the file and line of the error location and an optional message.
#define warning(format, ...) {\
        fprintf(stderr, "sl2cfoam: WARNING at file %s, line %d. Continuing...\n", __FILE__, __LINE__);\
        fprintf(stderr, "        : ");\
        fprintf(stderr, format, ##__VA_ARGS__);\
        fprintf(stderr, "\n");\
        }

// Raises an error when a GSL command failed.
#define check_gsl(status) {\
        if (status) {\
            error("GSL error code %d", status);\
        }\
        }

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_ERROR_H__*/