#include <stdio.h>
#include "pll.h"

// declarations
pll_rtree_t * pll_rtree_parse_newick(const char * filename);
void pll_rtree_show_ascii(const pll_rnode_t * root, int options);
char * pll_rtree_export_newick(const pll_rnode_t * root,
                               char * (*cb_serialize)(const pll_rnode_t *));
struct pll_rtree_s;
struct pll_rtree_t;

int cb_rfull_traversal(pll_rnode_t * node)
{
  return 1;
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

    pll_rtree_show_ascii(tree->root,
                         PLL_UTREE_SHOW_LABEL |
                         PLL_UTREE_SHOW_BRANCH_LENGTH |
                         PLL_UTREE_SHOW_CLV_INDEX);
    char * newick = pll_rtree_export_newick(tree->root,NULL);
    printf("%s\n", newick);
    free(newick);

    pll_operation_t * operations;
    pll_rnode_t ** travbuffer;
    double * branch_lengths;
    unsigned int * matrix_indices;
    unsigned int matrix_count, ops_count;

    travbuffer = (pll_rnode_t **)malloc(nodes_count * sizeof(pll_rnode_t *));

    branch_lengths = (double *)malloc(branch_count * sizeof(double));
    matrix_indices = (unsigned int *)malloc(branch_count * sizeof(int));
    operations = (pll_operation_t *)malloc(inner_nodes_count *
                                                    sizeof(pll_operation_t));
    unsigned int traversal_size;

    if (!pll_rtree_traverse(tree->root,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            cb_rfull_traversal,
                            travbuffer,
                            &traversal_size))
        printf("Error running Function pll_rtree_traverse() with root node as parameter");

    pll_rtree_create_operations(travbuffer,
                                traversal_size,
                                branch_lengths,
                                matrix_indices,
                                operations,
                                &matrix_count,
                                &ops_count);
    printf("Total operations: %d\n", ops_count);

    for (int i = 0; i < ops_count; ++i) {
        printf("Operation %d:\t", i);
        printf("  Parent CLV Index: %u\t", operations[i].parent_clv_index);
        printf("  Child1 CLV Index: %u\t", operations[i].child1_clv_index);
        printf("  Child2 CLV Index: %u\n", operations[i].child2_clv_index);
        printf("\n");
    }

    return 0;
}
