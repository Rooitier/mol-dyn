/**
 * @file    linkedlist.h
 * @brief   Generic doubly-linked list with insertion at beginning and removal.
 *
 * TODO
 *
 * @author  D. Trinchero (dario.trinchero@gmail.com)
 * @date    2019-09-19
 */

#ifndef LINKEDLIST_H
#define INKEDLIST_H

#include <stdlib.h>

/** the container structure for a linked list node */
typedef struct listnode ListNode;

/** a linked list node container */
struct listnode {
    void *data;
    ListNode *next;
    ListNode *prev;
};

/** the container structure for a linked list */
typedef struct linkedlist LinkedList;

/** a linked list container */
struct linkedlist {
    ListNode *head;
};

/**
 * Creates and returns a new (empty) list.
 *
 * @return      a pointer to the new list which was created
 */
LinkedList *lst_new_list(void);

/**
 * Creates and returns a new list node with the given data.
 *
 * @param[in]   data
 *     a pointer to the data which the new node is to contain
 * @return      a pointer to the new node which was created
 */
ListNode *lst_new_node(void *data);

/**
 * Links the given node at the beginning of the given list.
 *
 * @param[in]   list
 *     a pointer to the list to which to append the node
 * @param[in]   node
 *     a pointer to the node which is to be appended
 */
void lst_link(LinkedList *list, ListNode *node);

/**
 * Unlinks and returns the given node from the given list and returns it.
 *
 * @param[in]   list
 *     a pointer to the list from which to remove the node
 * @param[in]   node
 *     a pointer to the node which is to be removed
 * @return      a pointer to the unlinked node
 */
ListNode *lst_unlink(LinkedList *list, ListNode *node);

/**
 * Moves the given node from one list to another.
 *
 * @param[in]   old_lst
 *     a pointer to the list from which to remove the node
 * @param[in]   new_lst
 *     a pointer to the list to which to add the node
 * @param[in]   node
 *     a pointer to the node which is to be moved
 */
void lst_move(LinkedList *old_lst, LinkedList *new_lst, ListNode *node);

/**
 * Frees all of the nodes in the linked list, along with the data they store.
 *
 * @param[in]   list
 *     a pointer to the list which is to be freed
 * @param[in]   free_data
 *     a pointer to a function which dellocates the node data
 */
void lst_free_list(LinkedList *list, void (*free_data)(void *k));

/**
 * Frees a single list node, along with the data it stores.
 *
 * @param[in]   node
 *     a pointer to the node which is to be freed
 * @param[in]   free_data
 *     a pointer to a function which dellocates the node data
 */
void lst_free_node(ListNode *node, void (*free_data)(void *k));

#endif
