/**
 * @file    linkedlist.c
 * @brief   Implementations for generic doubly-linked list API.
 *
 * TODO
 *
 * @author  D. Trinchero (dario.trinchero@gmail.com)
 * @date    2019-09-19
 */

#include "error.h"
#include "linkedlist.h"

/* --- linked list interface ------------------------------------------------ */

LinkedList *lst_new_list(void)
{
    LinkedList *new = emalloc(sizeof(LinkedList));
    new->head = NULL;
    return new;
}

ListNode *lst_new_node(void *data)
{
    ListNode *new = emalloc(sizeof(ListNode));
    new->prev = NULL;
    new->next = NULL;
    new->data = data;
    return new;
}

void lst_link(LinkedList *list, ListNode *node)
{
    node->prev = NULL;
    node->next = list->head;
    if (list->head) list->head->prev = node;
    list->head = node;
}

ListNode *lst_unlink(LinkedList *list, ListNode *node)
{
    if (node->next) node->next->prev = node->prev;
    if (node->prev) node->prev->next = node->next;
    else list->head = node->next;
    node->prev = NULL;
    node->next = NULL;
    return node;
}

void lst_move(LinkedList *old_lst, LinkedList *new_lst, ListNode *node)
{
	lst_link(new_lst, lst_unlink(old_lst, node));
}

void lst_free_list(LinkedList *list, void (*free_data)(void *k))
{
    ListNode *p, *q;
    for (p = list->head; p; p = q) {
        q = p->next;
        lst_free_node(p, free_data);
    }
    free(list);
}

void lst_free_node(ListNode *node, void (*free_data)(void *k))
{
    (*free_data)(node->data);
    free(node);
}
