#ifndef JSW_AVLTREE_H
#define JSW_AVLTREE_H

/*
  AVL balanced tree library

    > Created (Julienne Walker): June 17, 2003
    > Modified (Julienne Walker): September 24, 2005

  This code is in the public domain. Anyone may
  use it or change it in any way that they see
  fit. The author assumes no responsibility for 
  damages incurred through use of the original
  code or any variations thereof.

  It is requested, but not required, that due
  credit is given to the original author and
  anyone who has modified the code through
  a header comment, such as this one.
*/
#include <stddef.h>

/* Opaque types */
typedef struct jsw_avltree jsw_avltree_t;
typedef struct jsw_avltrav jsw_avltrav_t;

/* User-defined item handling */
typedef int   (*cmp_f) ( const void *p1, const void *p2 );
typedef void *(*dup_f) ( void *p );
typedef void  (*rel_f) ( void *p );

int int_cmp_f ( const void *p1, const void *p2 );
void *int_dup_f( void *p );
void int_rel_f( void *p );

int uint_cmp_f ( const void *p1, const void *p2 );
void *uint_dup_f( void *p );
void uint_rel_f( void *p );

int uint64_cmp_f ( const void *p1, const void *p2 );
void *uint64_dup_f( void *p );
void uint64_rel_f( void *p );

/* AVL tree functions */
jsw_avltree_t *jsw_avlnew ( cmp_f cmp, dup_f dup, rel_f rel );
void           jsw_avldelete ( jsw_avltree_t *tree );
void          *jsw_avlfind ( jsw_avltree_t *tree, void *data );
int            jsw_avlinsert ( jsw_avltree_t *tree, void *data );
int            jsw_avlerase ( jsw_avltree_t *tree, void *data );
size_t         jsw_avlsize ( jsw_avltree_t *tree );

/* Traversal functions */
jsw_avltrav_t *jsw_avltnew ( void );
void           jsw_avltdelete ( jsw_avltrav_t *trav );
void          *jsw_avltfirst ( jsw_avltrav_t *trav, jsw_avltree_t *tree );
void          *jsw_avltlast ( jsw_avltrav_t *trav, jsw_avltree_t *tree );
void          *jsw_avltnext ( jsw_avltrav_t *trav );
void          *jsw_avltprev ( jsw_avltrav_t *trav );

#endif
