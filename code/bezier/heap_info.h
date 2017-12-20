#ifndef HEAP_INFO_H
#define HEAP_INFO_H
/*** heap_info.h - Defines the universal heap_info_t data structure which each
 *** heap can provide, and algorithms can use.  An algorithm which uses the
 *** universal definition of a heap data structure can use different heaps
 *** interchangeably for testing or comparison purposes.
 ***/
/*
 *  Shane Saunders
 */

/* Structure to be provided by heaps and used by algorithms. */
typedef struct heap_info {
    long (*delete_min)(void *heap);
    void (*insert)(void *heap, long node, float key);
    void (*decrease_key)(void *heap, long node, float newkey);
    long (*n)(void *heap);
    long (*key_comps)(void *heap);
    void *(*alloc)(long n_items);
    void (*free)(void *heap);
    void (*dump)(void *heap);
} heap_info_t;

#endif
