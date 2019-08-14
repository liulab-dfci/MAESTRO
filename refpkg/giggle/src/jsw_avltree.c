/*
  AVL balanced tree library

    > Created (Julienne Walker): June 17, 2003
    > Modified (Julienne Walker): September 24, 2005
*/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <err.h>

#include "jsw_avltree.h"

#ifndef HEIGHT_LIMIT
#define HEIGHT_LIMIT 64 /* Tallest allowable tree */
#endif

typedef struct jsw_avlnode {
  int                 balance; /* Balance factor */
  void               *data;    /* User-defined content */
  struct jsw_avlnode *link[2]; /* Left (0) and right (1) links */
} jsw_avlnode_t;

struct jsw_avltree {
  jsw_avlnode_t *root; /* Top of the tree */
  cmp_f          cmp;    /* Compare two items */
  dup_f          dup;    /* Clone an item (user-defined) */
  rel_f          rel;    /* Destroy an item (user-defined) */
  size_t         size;   /* Number of items (user-defined) */
};

struct jsw_avltrav {
  jsw_avltree_t *tree;               /* Paired tree */
  jsw_avlnode_t *it;                 /* Current node */
  jsw_avlnode_t *path[HEIGHT_LIMIT]; /* Traversal path */
  size_t         top;                /* Top of stack */
};

/* Two way single rotation */
#define jsw_single(root,dir) do {         \
  jsw_avlnode_t *save = root->link[!dir]; \
  root->link[!dir] = save->link[dir];     \
  save->link[dir] = root;                 \
  root = save;                            \
} while (0)

/* Two way double rotation */
#define jsw_double(root,dir) do {                    \
  jsw_avlnode_t *save = root->link[!dir]->link[dir]; \
  root->link[!dir]->link[dir] = save->link[!dir];    \
  save->link[!dir] = root->link[!dir];               \
  root->link[!dir] = save;                           \
  save = root->link[!dir];                           \
  root->link[!dir] = save->link[dir];                \
  save->link[dir] = root;                            \
  root = save;                                       \
} while (0)

/* Adjust balance before double rotation */
#define jsw_adjust_balance(root,dir,bal) do { \
  jsw_avlnode_t *n = root->link[dir];         \
  jsw_avlnode_t *nn = n->link[!dir];          \
  if ( nn->balance == 0 )                     \
    root->balance = n->balance = 0;           \
  else if ( nn->balance == bal ) {            \
    root->balance = -bal;                     \
    n->balance = 0;                           \
  }                                           \
  else { /* nn->balance == -bal */            \
    root->balance = 0;                        \
    n->balance = bal;                         \
  }                                           \
  nn->balance = 0;                            \
} while (0)

/* Rebalance after insertion */
#define jsw_insert_balance(root,dir) do {  \
  jsw_avlnode_t *n = root->link[dir];      \
  int bal = dir == 0 ? -1 : +1;            \
  if ( n->balance == bal ) {               \
    root->balance = n->balance = 0;        \
    jsw_single ( root, !dir );             \
  }                                        \
  else { /* n->balance == -bal */          \
    jsw_adjust_balance ( root, dir, bal ); \
    jsw_double ( root, !dir );             \
  }                                        \
} while (0)

/* Rebalance after deletion */
#define jsw_remove_balance(root,dir,done) do { \
  jsw_avlnode_t *n = root->link[!dir];         \
  int bal = dir == 0 ? -1 : +1;                \
  if ( n->balance == -bal ) {                  \
    root->balance = n->balance = 0;            \
    jsw_single ( root, dir );                  \
  }                                            \
  else if ( n->balance == bal ) {              \
    jsw_adjust_balance ( root, !dir, -bal );   \
    jsw_double ( root, dir );                  \
  }                                            \
  else { /* n->balance == 0 */                 \
    root->balance = -bal;                      \
    n->balance = bal;                          \
    jsw_single ( root, dir );                  \
    done = 1;                                  \
  }                                            \
} while (0)

static jsw_avlnode_t *new_node ( jsw_avltree_t *tree, void *data )
{
  jsw_avlnode_t *rn = (jsw_avlnode_t *)malloc ( sizeof *rn );

  if ( rn == NULL )
    err(1, "malloc error at new_node()");

  rn->balance = 0;
  rn->data = tree->dup ( data );
  rn->link[0] = rn->link[1] = NULL;

  return rn;
}

jsw_avltree_t *jsw_avlnew ( cmp_f cmp, dup_f dup, rel_f rel )
{
  jsw_avltree_t *rt = (jsw_avltree_t *)malloc ( sizeof *rt );

  if ( rt == NULL )
      err(1, "malloc error at jsw_avlnew()");

  rt->root = NULL;
  rt->cmp = cmp;
  rt->dup = dup;
  rt->rel = rel;
  rt->size = 0;

  return rt;
}

void jsw_avldelete ( jsw_avltree_t *tree )
{
  jsw_avlnode_t *it = tree->root;
  jsw_avlnode_t *save;

  /* Destruction by rotation */
  while ( it != NULL ) {
    if ( it->link[0] == NULL ) {
      /* Remove node */
      save = it->link[1];
      tree->rel ( it->data );
      free ( it );
    }
    else {
      /* Rotate right */
      save = it->link[0];
      it->link[0] = save->link[1];
      save->link[1] = it;
    }

    it = save;
  }

  free ( tree );
}

void *jsw_avlfind ( jsw_avltree_t *tree, void *data )
{
  jsw_avlnode_t *it = tree->root;

  while ( it != NULL ) {
    int cmp = tree->cmp ( it->data, data );

    if ( cmp == 0 )
      break;

    it = it->link[cmp < 0];
  }

  return it == NULL ? NULL : it->data;
}

int jsw_avlinsert ( jsw_avltree_t *tree, void *data )
{
  /* Empty tree case */
  if ( tree->root == NULL ) {
    tree->root = new_node ( tree, data );
    if ( tree->root == NULL )
      return 0;
  }
  else {
    jsw_avlnode_t head = {0}; /* Temporary tree root */
    jsw_avlnode_t *s, *t;     /* Place to rebalance and parent */
    jsw_avlnode_t *p, *q;     /* Iterator and save pointer */
    int dir;

    /* Set up false root to ease maintenance */
    t = &head;
    t->link[1] = tree->root;

    /* Search down the tree, saving rebalance points */
    for ( s = p = t->link[1]; ; p = q ) {
      dir = tree->cmp ( p->data, data ) < 0;
      q = p->link[dir];

      if ( q == NULL )
        break;
      
      if ( q->balance != 0 ) {
        t = p;
        s = q;
      }
    }

    p->link[dir] = q = new_node ( tree, data );
    if ( q == NULL )
      return 0;

    /* Update balance factors */
    for ( p = s; p != q; p = p->link[dir] ) {
      dir = tree->cmp ( p->data, data ) < 0;
      p->balance += dir == 0 ? -1 : +1;
    }

    q = s; /* Save rebalance point for parent fix */

    /* Rebalance if necessary */
    if ( abs ( s->balance ) > 1 ) {
      dir = tree->cmp ( s->data, data ) < 0;
      jsw_insert_balance ( s, dir );
    }

    /* Fix parent */
    if ( q == head.link[1] )
      tree->root = s;
    else
      t->link[q == t->link[1]] = s;
  }

  ++tree->size;

  return 1;
}

int jsw_avlerase ( jsw_avltree_t *tree, void *data )
{
  if ( tree->root != NULL ) {
    jsw_avlnode_t *it, *up[HEIGHT_LIMIT];
    int upd[HEIGHT_LIMIT], top = 0;
    int done = 0;

    it = tree->root;

    /* Search down tree and save path */
    for ( ; ; ) {
      if ( it == NULL )
        return 0;
      else if ( tree->cmp ( it->data, data ) == 0 )
        break;

      /* Push direction and node onto stack */
      upd[top] = tree->cmp ( it->data, data ) < 0;
      up[top++] = it;

      it = it->link[upd[top - 1]];
    }

    /* Remove the node */
    if ( it->link[0] == NULL || it->link[1] == NULL ) {
      /* Which child is not null? */
      int dir = it->link[0] == NULL;

      /* Fix parent */
      if ( top != 0 )
        up[top - 1]->link[upd[top - 1]] = it->link[dir];
      else
        tree->root = it->link[dir];

      tree->rel ( it->data );
      free ( it );
    }
    else {
      /* Find the inorder successor */
      jsw_avlnode_t *heir = it->link[1];
      void *save;
      
      /* Save this path too */
      upd[top] = 1;
      up[top++] = it;

      while ( heir->link[0] != NULL ) {
        upd[top] = 0;
        up[top++] = heir;
        heir = heir->link[0];
      }

      /* Swap data */
      save = it->data;
      it->data = heir->data;
      heir->data = save;

      /* Unlink successor and fix parent */
      up[top - 1]->link[up[top - 1] == it] = heir->link[1];

      tree->rel ( heir->data );
      free ( heir );
    }

    /* Walk back up the search path */
    while ( --top >= 0 && !done ) {
      /* Update balance factors */
      up[top]->balance += upd[top] != 0 ? -1 : +1;

      /* Terminate or rebalance as necessary */
      if ( abs ( up[top]->balance ) == 1 )
        break;
      else if ( abs ( up[top]->balance ) > 1 ) {
        jsw_remove_balance ( up[top], upd[top], done );

        /* Fix parent */
        if ( top != 0 )
          up[top - 1]->link[upd[top - 1]] = up[top];
        else
          tree->root = up[0];
      }
    }

    --tree->size;
  }

  return 1;
}

size_t jsw_avlsize ( jsw_avltree_t *tree )
{
  return tree->size;
}

jsw_avltrav_t *jsw_avltnew ( void )
{
  return malloc ( sizeof ( jsw_avltrav_t ) );
}

void jsw_avltdelete ( jsw_avltrav_t *trav )
{
  free ( trav );
}

/*
  First step in traversal,
  handles min and max
*/
static void *start ( jsw_avltrav_t *trav, jsw_avltree_t *tree, int dir )
{
  trav->tree = tree;
  trav->it = tree->root;
  trav->top = 0;

  /* Build a path to work with */
  if ( trav->it != NULL ) {
    while ( trav->it->link[dir] != NULL ) {
      trav->path[trav->top++] = trav->it;
      trav->it = trav->it->link[dir];
    }
  }

  return trav->it == NULL ? NULL : trav->it->data;
}

/*
  Subsequent traversal steps,
  handles ascending and descending
*/
static void *move ( jsw_avltrav_t *trav, int dir )
{
  if ( trav->it->link[dir] != NULL ) {
    /* Continue down this branch */
    trav->path[trav->top++] = trav->it;
    trav->it = trav->it->link[dir];

    while ( trav->it->link[!dir] != NULL ) {
      trav->path[trav->top++] = trav->it;
      trav->it = trav->it->link[!dir];
    }
  }
  else {
    /* Move to the next branch */
    jsw_avlnode_t *last;

    do {
      if ( trav->top == 0 ) {
        trav->it = NULL;
        break;
      }

      last = trav->it;
      trav->it = trav->path[--trav->top];
    } while ( last == trav->it->link[dir] );
  }

  return trav->it == NULL ? NULL : trav->it->data;
}

void *jsw_avltfirst ( jsw_avltrav_t *trav, jsw_avltree_t *tree )
{
  return start ( trav, tree, 0 ); /* Min value */
}

void *jsw_avltlast ( jsw_avltrav_t *trav, jsw_avltree_t *tree )
{
  return start ( trav, tree, 1 ); /* Max value */
}

void *jsw_avltnext ( jsw_avltrav_t *trav )
{
  return move ( trav, 1 ); /* Toward larger items */
}

void *jsw_avltprev ( jsw_avltrav_t *trav )
{
  return move ( trav, 0 ); /* Toward smaller items */
}

int int_cmp_f ( const void *p1, const void *p2 )
{
    int *a = (int *)p1;
    int *b = (int *)p2;

    return *a - *b;
}

void *int_dup_f( void *p )
{
    int *o = (int *)p;
    int *a = (int *)malloc(sizeof(int));
    if (a == NULL)
        err(1, "malloc error in int_dup_f");

    *a = *o;

    return a;
}

void int_rel_f( void *p )
{
    free(p);
}

int uint_cmp_f ( const void *p1, const void *p2 )
{
    uint32_t *a = (uint32_t *)p1;
    uint32_t *b = (uint32_t *)p2;

    if (*a < *b)
        return -1;
    else if (*a > *b)
        return 1;
    else 
        return 0;
}

void *uint_dup_f( void *p )
{
    uint32_t *o = (uint32_t *)p;
    uint32_t *a = (uint32_t *)malloc(sizeof(uint32_t));
    if (a == NULL)
        err(1, "malloc error in uint_dup_f()");

    *a = *o;

    return a;
}

void uint_rel_f( void *p )
{
    free(p);
}


int uint64_cmp_f ( const void *p1, const void *p2 )
{
    uint64_t *a = (uint64_t *)p1;
    uint64_t *b = (uint64_t *)p2;

    if (*a < *b)
        return -1;
    else if (*a > *b)
        return 1;
    else 
        return 0;
}

void *uint64_dup_f( void *p )
{
    uint64_t *o = (uint64_t *)p;
    uint64_t *a = (uint64_t *)malloc(sizeof(uint64_t));
    if (a == NULL)
        err(1, "malloc error in uint_dup_f()");

    *a = *o;

    return a;
}

void uint64_rel_f( void *p )
{
    free(p);
}


