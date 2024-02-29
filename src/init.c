#include <stdio.h>
#include "pll.h"
#include "init.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include <assert.h>

#define max(a, b) ((a) > (b) ? (a) : (b))

// declarations
pll_rtree_t * pll_rtree_parse_newick(const char * filename);
void pll_rtree_show_ascii(const pll_rnode_t * root, int options);
char * pll_rtree_export_newick(const pll_rnode_t * root,
                               char * (*cb_serialize)(const pll_rnode_t *));
void pll_rtree_destroy(pll_rtree_t * root,
                       void (*cb_destroy)(void *));
pll_partition_t * pll_partition_create(unsigned int tips,
                                       unsigned int clv_buffers,
                                       unsigned int states,
                                       unsigned int sites,
                                       unsigned int rate_matrices,
                                       unsigned int prob_matrices,
                                       unsigned int rate_cats,
                                       unsigned int scale_buffers,
                                       unsigned int attributes);
void pll_partition_destroy(pll_partition_t * partition);
int pll_set_tip_states(pll_partition_t * partition,
                       unsigned int tip_index,
                       const pll_state_t * map,
                       const char * sequence);

struct pll_rtree_s;
struct pll_rtree_t;
struct pll_partition_t;
struct pll_state_t;

const size_t N_CELLS = 1000;
// update this constant with a more realistic value
const size_t N_SITES = 10;
const int RATE_CATS = 1;
const size_t N_ALLELES = 3000;

char *xstrdup(const char *s) {
    if (s == NULL) {
        return NULL;
    }
    size_t len = strlen(s) + 1;
    char *dup = malloc(len);
    if (dup == NULL) {
        return NULL;
    }
    strcpy(dup, s);
    return dup;
}

void insert(ENTRY** tbl, ENTRY* entry) {
    hsearch(*entry, ENTER);
}

ENTRY* search(ENTRY** tbl, const char* key) {
    ENTRY query;
    query.key = (char*)key;

    ENTRY* res = hsearch(query, FIND);
    if (res != NULL) {
        return res;
    }
    return NULL;
}

void destroy(ENTRY** tbl) {
    hdestroy();
}

void print_binary(pll_state_t num) {
    for (size_t i = sizeof(unsigned long long) * 8 - 1; i >= 0; i--) {
        unsigned long long bit = (num >> i) & 1;
        printf("%llu", bit);
        if (i % 4 == 0) {
            printf(" ");
        }
    }
    printf("\n");
}

void print_sites(qfml_site_t ** site_data, size_t n_row) {
    for (size_t i = 0; i < n_row; ++i) {
        qfml_site_t site = (*site_data)[i];
        printf("%s\t%d", site.id, site.num_alleles);
        qfml_allele_t * alleles= site.alleles;
        for (size_t j = 0; j < site.num_alleles; ++j) {
            qfml_allele_t allele = alleles[j];
            printf("\tR-%d\t%e", allele.id, allele.prob);
        }
        printf("\n");
    }

}

// cmat operations
void init_cmat(unsigned int*** cmat, size_t n_row, size_t n_col) {
    *cmat = (unsigned int**)malloc(n_row * sizeof(unsigned int*));
    (*cmat)[0] = (unsigned int*)malloc(n_row * n_col * sizeof(unsigned int));

    for (size_t i = 1; i < n_row; ++i) {
        (*cmat)[i] = (*cmat)[0] + i * n_col;
    }
}

void resize_cmat_horiz(unsigned int*** cmat, size_t n_row, size_t* n_col) {
    size_t n_col_new = *n_col * 2;
    *cmat = (unsigned int**)realloc(*cmat, n_row * n_col_new * sizeof(unsigned int));
    if (*cmat != NULL) {
        *n_col = n_col_new;
    } else {
        perror("(resize_cmat_horiz) failed resize operation");
        exit(1);
    }
}

void resize_cmat_vert(unsigned int*** cmat, size_t* n_row, size_t n_col) {
    size_t n_row_new = *n_row * 2;
    *cmat = (unsigned int**)realloc(*cmat, n_row_new * sizeof(unsigned int*));
    if (*cmat != NULL) {
        for (size_t i = *n_row; i < n_row_new; ++i) {
            (*cmat)[i] = (unsigned int*)malloc(n_col * sizeof(unsigned int));
        }
        *n_row = n_row_new;
    } else {
        perror("(resize_cmat_vert) failed resize operation");
        exit(1);
    }
}

void free_cmat(unsigned int*** cmat) {
    free((*cmat)[0]);
    free(*cmat);
}

void print_cmat(unsigned int*** cmat, size_t n_row, size_t n_col, char*** cell_ids) {
    for (size_t i = 0; i < n_row; ++i) {
        printf("%s\t", (*cell_ids)[i]);
        for (size_t j = 0; j < n_col; ++j) {
            printf("%u\t", (*cmat)[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char * argv[]) {
    printf("Parsing newick formatted tree...\n");
    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
    pll_partition_t * partition;

    pll_rtree_t * tree = pll_rtree_parse_newick(argv[1]);
    if (!tree) {
        printf("Error parsing the tree\n");
        exit(1);
    }

    tip_nodes_count = tree->tip_count;

    /* compute and show node count information */
    inner_nodes_count = tip_nodes_count - 1;
    nodes_count = inner_nodes_count + tip_nodes_count;
    branch_count = nodes_count - 1;

    printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
    printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
    printf("Total number of nodes in tree: %d\n", nodes_count);
    printf("Number of branches in tree: %d\n", branch_count);

    // store <tip label, tip clv index> kv pairs in a libc hash table
    ENTRY* ht1 = NULL;
    hcreate(tip_nodes_count);
    unsigned int * data = (unsigned int *)malloc(tip_nodes_count * sizeof(unsigned int));
    for (size_t i = 0; i < tip_nodes_count; ++i) {
        data[i] = tree->nodes[i]->clv_index;
        ENTRY entry;
#ifdef __APPLE__
        entry.key = xstrdup(tree->nodes[i]->label);
#else
        entry.key = tree->nodes[i]->label;
#endif
        entry.data = (void *)(data+i);
        insert(&ht1, &entry);
    }

    // parse in chr_mat.tsv
    FILE* cmat_file = fopen(argv[2], "r");

    if (cmat_file == NULL) {
        printf("Error opening chr mat file.\n");
        exit(1);
    }

    int first_ln = 1;
    char* ln = NULL;
    size_t ln_size = 0;

    unsigned int** cmat = NULL;
    init_cmat(&cmat, N_CELLS, N_SITES);
    size_t n_row = N_CELLS;
    size_t n_col = N_SITES;

    char** cell_ids = (char**)malloc(n_row * sizeof(char*));
    int cell_idx = 0;
    int site_idx = 0;

    while (getline(&ln, &ln_size, cmat_file) != -1) {
        // skip the first line (i.e., header)
        if (first_ln) {
            first_ln = 0;
            continue;
        }

        // store tokens first
        char** tokens = (char**)malloc(n_col * sizeof(char*));
        int idx = 0;
        char* tok = strtok(ln, "\t\n");
        while (tok != NULL) {
            // re-size if needed
            // after iterating through the first line, no more vertical resizing will be executed.
            if (idx == n_col) {
                // TODO: change to a more informative message
                printf("expanding horizontally\n");
                resize_cmat_horiz(&cmat, n_row, &n_col);
                tokens = (char**)realloc(tokens, n_col * sizeof(char*));
            }
            tokens[idx] = strdup(tok);
            idx++;
            tok = strtok(NULL, "\t\n");
        }
        free(tok);
        //printf("Loading Cell ID: %s\n", tokens[0]);

        if (cell_idx == n_row) {
            // cell_ids is re-used
            // TODO: change to a more informative message
            printf("expanding vertically\n");
            resize_cmat_vert(&cmat, &n_row, n_col);
            cell_ids = (char**)realloc(cell_ids, n_row * sizeof(char*));
        }

        cell_ids[cell_idx] = strdup(tokens[0]);
        for (size_t i = 1; i < idx; ++i) {
            unsigned int allele_id = atoi(tokens[i] + 2);
            cmat[cell_idx][i - 1] = allele_id;
            //printf("Original allele info: %s\n", tokens[i]);
            //printf("Stored allele value: %d\n", cmat[cell_idx][i - 1]);
        }
        cell_idx++;
        site_idx = max(site_idx, idx - 1);
        free(tokens);
    }
    // store new cell / site numbers
    size_t num_cells = cell_idx;
    size_t num_sites = site_idx;
    //print_cmat(&cmat, num_cells, num_sites, &cell_ids);

    FILE* allele_prob_file = fopen(argv[3], "r");
    if (allele_prob_file == NULL) {
        printf("Error opening chr mat file.\n");
        exit(1);
    }

    first_ln = 1;
    ln = NULL;
    ln_size = 0;

    qfml_site_t * site_data = (qfml_site_t *)malloc(num_sites *
                                                  sizeof(qfml_site_t));
    site_idx = 0;
    char * prev_site = NULL;
    size_t alleles_size = N_ALLELES;
    qfml_allele_t * alleles = (qfml_allele_t *)malloc(alleles_size * sizeof(qfml_allele_t));
    int allele_idx = 0;

    while (getline(&ln, &ln_size, allele_prob_file) != -1) {
        if (first_ln) {
            first_ln = 0;
            continue;
        }
        int idx = 0;
        char ** tokens = (char**)malloc(3 * sizeof(char*));
        char* tok = strtok(ln, "\t\n");
        while (tok != NULL) {
            tokens[idx] = strdup(tok);
            idx++;
            tok = strtok(NULL, "\t\n");
        }
        free(tok);
        if (prev_site == NULL) {
            prev_site = strdup(tokens[0]);
        } else {
            if (strcmp(prev_site, tokens[0]) != 0) {
                qfml_site_t site;
                site.num_alleles = allele_idx;
                site.alleles = alleles;
                site.id = prev_site;
                site_data[site_idx] = site;
                site_idx++;
                prev_site = strdup(tokens[0]);
                alleles_size = N_ALLELES;
                free(alleles);
                alleles = (qfml_allele_t *)malloc(alleles_size * sizeof(qfml_allele_t));
                allele_idx = 0;
            }
        }
        // resize if needed
        if (allele_idx == alleles_size) {
//            printf("expanding alleles vector\n");
            alleles_size *= 2;
            alleles = (qfml_allele_t *)realloc(alleles, alleles_size * sizeof(qfml_allele_t));
        }
        qfml_allele_t allele;
        allele.id = atoi(tokens[1] + 2);
        allele.prob = strtod(tokens[2], NULL);
        alleles[allele_idx] = allele;
        allele_idx++;
        free(tokens);
    }
    // shallow copy
    qfml_site_t site;
    site.num_alleles = allele_idx;
    site.alleles = alleles;
    site.id = prev_site;
    site_data[site_idx] = site;

    // sanity check
    assert(site_idx == num_sites);

    //print_sites(&site_data, num_sites);

    ENTRY* ht2 = NULL;
    hcreate(num_sites);

    for (size_t i = 0; i < num_sites; ++i) {
        ENTRY entry;

#ifdef __APPLE__
        entry.key = xstrdup(site_data[i].id);
#else
        entry.key = strdup(site_data[i].id);
#endif
        if (entry.key == NULL) {
            perror("(xstrdup) failed to allocate memory");
            exit(1);
        }
        entry.data = (void *)(site_data+i);
        insert(&ht2, &entry);
    }

//    // test out ht2 search operation
//    printf("%s\t", site_data[0].id);
//    ENTRY * res = search(&ht2, site_data[0].id);
//    printf("%d\n", ((qfml_site_t *)res->data)->num_alleles);
//
//    // test out ht1 search operation
//    printf("%s\t", tree->nodes[0]->label);
//    res = search(&ht1, tree->nodes[0]->label);
//    printf("%d\n", *((unsigned int *)(res->data)));

    // create pll partition
    // experiment with a single site first
    partition = pll_partition_create(tip_nodes_count,
                                     inner_nodes_count,
                                     // number of states
                                     // should add the null state
                                     site_data[0].num_alleles + 1,
                                     snprintf(NULL, 0, "%d", site_data[0].num_alleles),
                                     1,
                                     branch_count,
                                     RATE_CATS,
                                     inner_nodes_count,
                                     PLL_ATTRIB_ARCH_CPU);

    // set tip chars
    assert(tip_nodes_count == num_cells);
    for (size_t i = 0; i < tip_nodes_count; ++i) {
        ENTRY * h1_res = search(&ht1, cell_ids[i]);
        if (!h1_res) {
            fprintf(stderr, "Fatal error: cell with ID %s doesn't appear in the tree.\n", cell_ids[i]);
            exit(1);
        }
        unsigned int tip_clv_index = *((unsigned int *)(h1_res->data));
        unsigned int cell_state = cmat[i][0];
        size_t len = snprintf(NULL, 0, "%d", cell_state);
        char *cell_state_s = (char *)malloc(len + 1);
        if (cell_state_s == NULL) {
            perror("(snprintf) unable to allocated memory\n");
            exit(1);
        }
        sprintf(cell_state_s, "%d", cell_state);
        // sanity check
//        for (size_t j = 0; cell_state_s[j] != '\0'; ++j) {
//            pll_state_t cell_pll_state = qfml_map_allele[cell_state_s[j]];
//            printf("%d\t", cell_state_s[j]);
//            print_binary(cell_pll_state);
//            printf("%llu\n", cell_pll_state);
//        }
        pll_set_tip_states(partition, tip_clv_index, qfml_map_allele, cell_state_s);
        free(cell_state_s);
    }
    free(data);

    // destroy all structures before exiting
    pll_partition_destroy(partition);
    free_cmat(&cmat);
    pll_rtree_destroy(tree,NULL);
    return 0;
}
