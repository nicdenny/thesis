/*** Fibonacci Heap Implementation ***/
/*
 *   Shane Saunders
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#if FHEAP_DUMP
#include <stdio.h>
#endif
#include "fheap.h"



/*** Prototypes of functions only visible within this file. ***/
void fh_dump_nodes(fheap_node_t *ptr, long level);
void fh_meld(fheap_t *h, fheap_node_t *tree_list);



/*** Definitions for functions that are visible outside this file. ***/

/* fh_alloc() - creates and and returns a pointer to a F-heap which can contain
 * up to max_nodes nodes.
 */
fheap_t *fh_alloc(long max_nodes)
{
    fheap_t *h;
#if FHEAP_DUMP
printf("alloc, ");
#endif
 
    /* Create the heap. */
    h = (fheap_t *)malloc(sizeof(fheap_t));
 
    h->max_trees = (long)(1.0 + 1.44 * log(max_nodes)/log(2.0));
    h->max_nodes = max_nodes;
    h->trees = (fheap_node_t **)calloc(h->max_trees, sizeof(fheap_node_t *));
    h->nodes = (fheap_node_t **)calloc(max_nodes, sizeof(fheap_node_t *));
    h->n = 0;

    /* The value of the heap helps to keep track of the maximum rank while
     * nodes are inserted or deleted.
     */
    h->value = 0;

    /* For experimental purposes, we keep a count of the number of key
     * comparisons.
     */
    h->key_comps = 0;

#if FHEAP_DUMP
printf("alloc-exited, ");
#endif
    return h; 
}


/* fh_free() - destroys the heap pointed to by h, freeing up any space that was
 * used by it.
 */
void fh_free(fheap_t *h)
{
    long i;
    
#if FHEAP_DUMP
printf("free, ");
#endif

    for(i = 0; i < h->max_nodes; i++) {
        free(h->nodes[i]);
    }

    free(h->nodes);
    free(h->trees);
    free(h);
    
#if FHEAP_DUMP
printf("free-exited, ");
#endif
}


/* fh_insert() - creates and inserts New a node representing vertex_no with key
 * k into the heap h.
 */
void fh_insert(fheap_t *h, long vertex_no, float k)
{
    fheap_node_t *New;

#if FHEAP_DUMP
printf("insert, ");
#endif

    /* Create an initialise the New node. */
    New = (fheap_node_t*)malloc(sizeof(fheap_node_t));
    New->child = NULL;
    New->left = New->right = New;
    New->rank = 0;
    New->vertex_no = vertex_no;
    New->key = k;

    /* Maintain a pointer vertex_no's New node in the heap. */
    h->nodes[vertex_no] = New;

    /* Meld the New node into the heap. */
    fh_meld(h, New);

    /* Update the heaps node count. */
    h->n++;

#if FHEAP_DUMP
printf("insert-exited, ");
#endif
}


/* fh_delete_min() - deletes the minimum node from the heap pointed to by h and
 * returns its vertex number.
 */
long fh_delete_min(fheap_t *h)
{
    fheap_node_t *min_node, *child, *next;
    float k, k2;
    long r, v, vertex_no;

#if FHEAP_DUMP
printf("delete_min, ");
#endif

    /* First we determine the maximum rank in the heap. */
    v = h->value;
    r = -1;
    while(v) {
        v = v >> 1;
        r++;
    };

    /* Now determine which root node is the minimum. */
    min_node = h->trees[r];
    k = min_node->key;
    while(r > 0) {
        r--;
        next = h->trees[r];
        if(next) {
            if((k2 = next->key) < k) {
                k = k2;
                min_node = next;
            }
            h->key_comps++;
        }
    }

    /* We remove the minimum node from the heap but keep a pointer to it. */
    r = min_node->rank;
    h->trees[r] = NULL;
    h->value -= (1 << r);

    child = min_node->child;
    if(child) fh_meld(h, child);

    /* Record the vertex no of the old minimum node before deleting it. */
    vertex_no = min_node->vertex_no;
    h->nodes[vertex_no] = NULL;
    free(min_node);
    h->n--;

#if FHEAP_DUMP
printf("delete_min-exited, ");
#endif

    return vertex_no;
}


/* fh_decrease_key() - decreases the key used for vertex, vertex_no, to
 * New_value.  No check is made to ensure that New_value is in-fact less than
 * the current value so it is up to the user of this function to ensure that
 * it is.
 */
void fh_decrease_key(fheap_t *h, long vertex_no, float New_value)
{
    fheap_node_t *cut_node, *parent, *New_roots, *r, *l;
    long prev_rank;

#if FHEAP_DUMP
printf("decrease_key on vn = %d, ", vertex_no);
#endif

    /* Obtain a pointer to the decreased node and its parent then decrease the
     * nodes key.
     */
    cut_node = h->nodes[vertex_no];
    parent = cut_node->parent;
    cut_node->key = New_value;

    /* No reinsertion occurs if the node changed was a root. */
    if(!parent) {
#if FHEAP_DUMP
printf("decrease_key-exited, ");
#endif
        return;
    }

    /* Update the left and right pointers of cut_node and its two neighbouring
     * nodes.
     */
    l = cut_node->left;
    r = cut_node->right;
    l->right = r;
    r->left = l;
    cut_node->left = cut_node->right = cut_node;

    /* Initially the list of New roots contains only one node. */
    New_roots = cut_node;

    /* While there is a parent node that is marked a cascading cut occurs. */
    while(parent && parent->marked) {

        /* Decrease the rank of cut_node's parent an update its child pointer.
         */
        parent->rank--;
        if(parent->rank) {
            if(parent->child == cut_node) parent->child = r;
        }
        else {
            parent->child = NULL;
        }

        /* Update the cut_node and parent pointers to the parent. */
        cut_node = parent;
        parent = cut_node->parent;

        /* Update the left and right pointers of cut_nodes two neighbouring
         * nodes.
         */
        l = cut_node->left;
        r = cut_node->right;
        l->right = r;
        r->left = l;

        /* Add cut_node to the list of nodes to be reinserted as New roots. */
        l = New_roots->left;
        New_roots->left = l->right = cut_node;
        cut_node->left = l;
        cut_node->right = New_roots;
        New_roots = cut_node;
    }

    /* If the root node is being relocated then update the trees[] array.
     * Otherwise mark the parent of the last node cut.
     */
    if(!parent) {
        prev_rank = cut_node->rank + 1;
        h->trees[prev_rank] = NULL;
        h->value -= (1 << prev_rank);
    }
    else {
        /* Decrease the rank of cut_node's parent an update its child pointer.
         */
        parent->rank--;
        if(parent->rank) {
            if(parent->child == cut_node) parent->child = r;
        }
        else {
            parent->child = NULL;
        }

        parent->marked = 1;
    }

    /* Meld the New roots into the heap. */
    fh_meld(h, New_roots);

#if FHEAP_DUMP
printf("decrease_key-exited, ");
#endif
}



/*** Definitions of functions that are only visible within this file. ***/

/* fh_meld() - melds  the linked list of trees pointed to by *tree_list into
 * the heap pointed to by h.
 */
void fh_meld(fheap_t *h, fheap_node_t *tree_list)
{
    fheap_node_t *first, *next, *node_ptr, *New_root, *temp, *temp2, *lc, *rc;
    long r;

#if FHEAP_DUMP
printf("meld: ");
#endif

    /* We meld each tree in the circularly linked list back into the root level
     * of the heap.  Each node in the linked list is the root node of a tree.
     * The circularly linked list uses the sibling pointers of nodes.  This
     *  makes melding of the child nodes from a delete_min operation simple.
     */
    node_ptr = first = tree_list;

    do {

#if FHEAP_DUMP
printf("%d, ", node_ptr->vertex_no);
#endif

        /* Keep a pointer to the next node and remove sibling and parent links
         * from the current node.  node_ptr points to the current node.
         */
        next = node_ptr->right;
        node_ptr->right = node_ptr->left = node_ptr;
        node_ptr->parent = NULL;

        /* We merge the current node, node_ptr, by inserting it into the
         * root level of the heap.
         */
        New_root = node_ptr;
        r = node_ptr->rank;

        /* This loop inserts the New root into the heap, possibly restructuring
         * the heap to ensure that only one tree for each degree exists.
         */
        do {

            /* Check if there is already a tree of degree r in the heap.
             * If there is then we need to link it with New_root so it will be
             * reinserted into a New place in the heap.
             */
            if((temp = h->trees[r])) {

	        /* temp will be linked to New_root and relocated so we no
                 * longer will have a tree of degree r.
                 */
                h->trees[r] = NULL;
                h->value -= (1 << r);

	        /* Swap temp and New_root if necessary so that New_root always
                 * points to the root node which has the smaller key of the
                 * two.
                 */
	        if(temp->key < New_root->key) {
                    temp2 = New_root;
                    New_root = temp;
                    temp = temp2;
                }
                h->key_comps++;

                /* Link temp with New_root, making sure that sibling pointers
                 * get updated if rank is greater than 0.  Also, increase r for
                 * the next pass through the loop since the rank of New has
                 * increased.
                 */
	        if(r++ > 0) {
                    rc = New_root->child;
                    lc = rc->left;
                    temp->left = lc;
                    temp->right = rc;
                    lc->right = rc->left = temp;
                }
                New_root->child = temp;
                New_root->rank = r;
                temp->parent = New_root;
                temp->marked = 0;
            }
            /* Otherwise if there is not a tree of degree r in the heap we
             * allow New_root, which possibly carries moved trees in the heap,
             * to be a tree of degree r in the heap.
             */
            else {

                h->trees[r] = New_root;
                h->value += (1 << r);;

                /* NOTE:  Because New_root is now a root we ensure it is
                 *        marked.
                 */
                New_root->marked = 1;
	    }

            /* Note that temp will be NULL if and only if there was not a tree
             * of degree r.
             */
        } while(temp);

        node_ptr = next;

    } while(node_ptr != first);

#if FHEAP_DUMP
printf("meld-exited, ");
#endif
}



/*** Debugging functions ***/

/* Recursively print the nodes of a Fibonacci heap. */
#if FHEAP_DUMP
void fh_dump_nodes(fheap_node_t *ptr, long level)
{
     fheap_node_t *child_ptr, *partner;
     long i, ch_count;

     /* Print leading whitespace for this level. */
     for(i = 0; i < level; i++) printf("   ");

     printf("%d(%ld)[%d]\n", ptr->vertex_no, ptr->key, ptr->rank);
     
     if((child_ptr = ptr->child)) {
	 child_ptr = ptr->child->right;
	 
         ch_count = 0;

         do {
             fh_dump_nodes(child_ptr, level+1);
	     if(child_ptr->dim > ptr->dim) {
                 for(i = 0; i < level+1; i++) printf("   ");
		 printf("error(dim)\n");  exit(1);
	     }
	     if(child_ptr->parent != ptr) {
                 for(i = 0; i < level+1; i++) printf("   ");
		 printf("error(parent)\n");
	     }
             child_ptr = child_ptr->right;
	     ch_count++;
         } while(child_ptr != ptr->child->right);

         if(ch_count != ptr->dim) {
	     for(i = 0; i < level; i++) printf("   ");
             printf("error(ch_count)\n");  exit(1);
         }
     }
     else { 
         if(ptr->dim != 0) {
             for(i = 0; i < level; i++) printf("   ");
	     printf("error(dim)\n"); exit(1);
	 }
     }

}
#endif

/* Print out a Fibonacci heap. */
#if FHEAP_DUMP
void fh_dump(fheap_t *h)
{
    long i;
    fheap_node_t *ptr;

    printf("\n");
    printf("value = %d\n", h->value);
    printf("array entries 0..max_trees =");
    for(i=0; i<h->max_trees; i++) {
        printf(" %d", h->trees[i] ? 1 : 0 );
    }
    printf("\n\n");
    for(i=0; i<h->max_trees; i++) {
        if((ptr = h->trees[i])) {
            printf("tree %d\n\n", i);
            fh_dump_nodes(ptr, 0);
	    printf("\n");
        }
    }
    fflush(stdout);
}
#endif


/*** Implement the universal heap structure type ***/

/* Fibonacci heap wrapper functions. */

long _fh_delete_min(void *h) {
    return fh_delete_min((fheap_t *)h);
}

void _fh_insert(void *h, long v, float k) {
    fh_insert((fheap_t *)h, v, k);
}

void _fh_decrease_key(void *h, long v, float k) {
    fh_decrease_key((fheap_t *)h, v, k);
}

long _fh_n(void *h) {
    return ((fheap_t *)h)->n;
}

long _fh_key_comps(void *h) {
    return ((fheap_t *)h)->key_comps;
}

void *_fh_alloc(long n) {
    return fh_alloc(n);
}

void _fh_free(void *h) {
    fh_free((fheap_t *)h);
}

void _fh_dump(void *h) {
#if FHEAP_DUMP
    fh_dump((fheap_t *)h);
#endif
}

/* Fibonacci heap info. */
const heap_info_t FHEAP_info = {
    _fh_delete_min,
    _fh_insert,
    _fh_decrease_key,
    _fh_n,
    _fh_key_comps,
    _fh_alloc,
    _fh_free,
    _fh_dump
};



/*
int main()
{
	fheap_t* H;
	
	H=fh_alloc(10000);
	fh_insert(H, 1, 123);
	fh_insert(H, 2, 1999.8);
	fh_insert(H, 3, 2000);
	fh_insert(H, 4, 200);
	long a=fh_delete_min(H);
	printf("%ld\n",a);
	
	fh_decrease_key(H, 3, 1999.9);

	a=fh_delete_min(H);
	printf("%ld %d\n",a,H->n);
	a=fh_delete_min(H);
	printf("%ld %d\n",a,H->n);
	a=fh_delete_min(H);
	printf("%ld %d\n",a,H->n);
	
}
*/


