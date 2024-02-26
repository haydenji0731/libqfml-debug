#include <stdio.h>
#include "pll.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// declarations
pll_rtree_t * pll_rtree_parse_newick(const char * filename);
void pll_rtree_show_ascii(const pll_rnode_t * root, int options);
char * pll_rtree_export_newick(const pll_rnode_t * root,
                               char * (*cb_serialize)(const pll_rnode_t *));
struct pll_rtree_s;
struct pll_rtree_t;

const size_t N_CELLS = 1000;
// update this constant with a more realistic value
const size_t N_SITES = 5;

char* convert_base(int n, int b) {
    if (n == 0) {
        char* base4_str = (char*)malloc(2 * sizeof(char));
        base4_str[0] = '0';
        base4_str[1] = '\0';
        return base4_str;
    }

    size_t mem_size = 3;
    char* out_str = (char*)malloc(mem_size * sizeof(char));
    size_t str_len = 0;
    int i = 0;

    while (n > 0) {
        // bound to add a character to the representation st
        int remainder = n % b;
        if (str_len + 1 == mem_size) {
            // resize & re-alloc memory
            mem_size *= 2;
            out_str = (char*)realloc(out_str, (mem_size) * sizeof(char));
        }

        out_str[i++] = remainder + '0';
        n /= b;
        str_len += 1;
    }
    out_str[i++] = '\0';
    return out_str;
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
        perror("failed resize operation.");
        exit(1);
    }
}

void resize_cmat_vert(unsigned int*** cmat, size_t* n_row, size_t n_col) {
    size_t n_row_new = *n_row * 2;
    *cmat = (unsigned int**)realloc(*cmat, n_row_new * n_col * sizeof(unsigned int));
    if (*cmat != NULL) {
        *n_row = n_row_new;
    } else {
        perror("failed resize operation.");
        exit(1);
    }
}

void free_cmat(unsigned int*** cmat) {
    free((*cmat)[0]);
    free(*cmat);
}

int main(int argc, char * argv[]) {
    printf("Parsing newick formatted tree...\n");
    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;

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

//    pll_rtree_show_ascii(tree->root,
//                         PLL_UTREE_SHOW_LABEL |
//                         PLL_UTREE_SHOW_BRANCH_LENGTH |
//                         PLL_UTREE_SHOW_CLV_INDEX);
//    char * newick = pll_rtree_export_newick(tree->root,NULL);
//    printf("%s\n", newick);
//    free(newick);

    // parse in chr_mat.tsv
    FILE* cmat_file = fopen(argv[2], "r");

    if (cmat_file == NULL) {
        printf("Error opening chr mat file.\n");
        exit(1);
    }

    int first_ln = 1;
    char* ln = NULL;
    size_t ln_size = 0;
    size_t cell_ids_size = N_CELLS;
    char** cell_ids = (char**)malloc(cell_ids_size * sizeof(char*));
    int cell_idx = 0;

    unsigned int** cmat = NULL;
    init_cmat(&cmat, N_CELLS, N_SITES);
    size_t n_row = N_CELLS;
    size_t n_col = N_SITES;

    while (getline(&ln, &ln_size, cmat_file) != -1) {
        // skip the first line (i.e., header)
        if (first_ln == 1) {
            first_ln = 0;
            continue;
        }
        // store tokens first
        size_t tokens_size = N_SITES;
        char** tokens = (char**)malloc(tokens_size * sizeof(char*));
        int idx = 0;
        char* tok = strtok(ln, "\t\n");
        while (tok != NULL) {
            // re-size if needed
            // after iterating through the first line, no more vertical resizing will be executed.
            if (idx == tokens_size) {
                // TODO: is this the best implementation?
                tokens_size *= 2;
                resize_cmat_horiz(&cmat, n_row, &n_col);
                tokens = (char**)realloc(tokens, tokens_size * sizeof(char*));
            }
            tokens[idx] = strdup(tok);
            idx++;
            tok = strtok(NULL, "\t\n");
        }
        printf("Loading Cell ID: %s\n", tokens[0]);
        if (cell_idx == cell_ids_size) {
            cell_ids_size *= 2;
            cell_ids = (char**)realloc(cell_ids, cell_ids_size * sizeof(char*));
            resize_cmat_vert(&cmat, &n_row, n_col);
        }
        cell_ids[cell_idx] = strdup(tokens[0]);
        for (int i = 1; i < idx; ++i) {
            unsigned int site = atoi(tokens[i] + 2);
            cmat[cell_idx][i] = site;
            printf("%s", tokens[i]);
            printf("\t%d\n", cmat[cell_idx][i]);
        }


    }

    free_cmat(&cmat);
    return 0;
}
