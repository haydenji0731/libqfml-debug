#include <stdio.h>
#include "pll.h"

// declarations
pll_rtree_t * pll_rtree_parse_newick(const char * filename);
void pll_rtree_show_ascii(const pll_rnode_t * root, int options);
char * pll_rtree_export_newick(const pll_rnode_t * root,
                               char * (*cb_serialize)(const pll_rnode_t *));
struct pll_rtree_s;
struct pll_rtree_t;

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

    return 0;
}
