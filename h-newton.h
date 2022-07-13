
#include "real.h"

typedef struct Parameters {
    real m;  // central mass
    real q_r, p_r, q_phi, p_phi;  // coordinates & momenta
    real h0;  // stored initial value of Hamiltonian
} parameters;
