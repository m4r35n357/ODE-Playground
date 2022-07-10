
#include "real.h"

typedef struct {
    real m;  // central mass
    real q_r, p_r, q_phi, p_phi;  // coordinates & momenta
    real h0;  // stored initial value of Hamiltonian
} parameters;

//void *get_p (int argc, char **argv, int va_begin);
