#ifndef INIT_H
#define INIT_H

#include "pll.h"

// pll_state_t = unsigned long (more bits than unsigned int from libpll1)
struct pll_state_t;

const pll_state_t qfml_map_allele[256];


typedef struct {
    unsigned int id;
    double prob;
} qfml_allele_t;

typedef struct {
    char * id;
    unsigned int num_alleles;
    qfml_allele_t * alleles;
} qfml_site_t;

#endif // INIT_H
