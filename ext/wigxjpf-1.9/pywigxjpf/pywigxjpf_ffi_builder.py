from cffi import FFI
ffibuilder = FFI()

# When compiled by 'python pywigxjpf_ffi_builder.py' according to
# Makefile (README), we are in this directory, so need a ../ prefix.
# When compiling from setup.py, the working directory is the parent.

updir = ""
if __name__ == "__main__": # "__cffi__" when compiling from setup
    updir = "../"

ffibuilder.set_source("pywigxjpf_ffi",
   r""" /* Compiled code.  All is in the library. */
        /* The header is included to check the shared prototypes. */
        #include "wigxjpf.h"
    """,
    include_dirs=[updir + 'inc/'],
    extra_objects=[updir + 'lib/libwigxjpf.a'])

ffibuilder.cdef("""
    /* Declarations shared from C to Python. */
    double wig3jj(int two_j1, int two_j2, int two_j3,
                  int two_m1, int two_m2, int two_m3);
    double wig6jj(int two_j1, int two_j2, int two_j3,
                  int two_m1, int two_m2, int two_m3);
    double wig9jj(int two_j1, int two_j2, int two_j3,
                  int two_j4, int two_j5, int two_j6,
                  int two_j7, int two_j8, int two_j9);
    void wig_table_init(int max_two_j, int wigner_type);
    void wig_table_free(void);
    void wig_temp_init(int max_two_j);
    void wig_temp_free(void);
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
