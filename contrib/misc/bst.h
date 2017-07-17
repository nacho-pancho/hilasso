/*
 * BST.h
 */
 
#include <cstdlib>

template <class Type>
struct bNode {
    Type key;
    bNode *link[2];
    bNode(){
        link[0] = NULL;
        link[1] = NULL;
    }
};
 
template <class Type>
class BST {
private:
  int count; 
    bNode<Type> *root;
 
    void clear(bNode<Type> *&ptr);
    bool search(const Type &item, bNode<Type> *&curr, bNode<Type> *&prev, bool &lr) const;
    inline Type inOrder(bNode<Type> *ptr) const;
    inline int subNodes(bNode<Type>* const &ptr) const;
    int height(bNode<Type>* const &ptr) const;
    bNode<Type>* minmax(bNode<Type> *ptr, const bool &lr) const;
public:
    BST();
    ~BST();
    void clear();
    bool isEmpty() const;
    bool insert(const Type &item);
    bool remove(const Type &item);
    bool search(const Type &item, Type *&ptr) const;
    Type min() const;
    Type max() const;
    int size() const;
    int height() const;
};
 
#include "BST.cpp"

 
template <class Type>
void BST<Type>::clear(bNode<Type> *&ptr)
{
    if (ptr != NULL) {
        clear(ptr->link[0]);
        clear(ptr->link[1]);
        delete ptr;
    }
}
 
template <class Type>
bool BST<Type>::search(const Type &item, bNode<Type> *&curr, bNode<Type> *&prev, bool &lr) const
{
    while (curr != NULL) {
        if (item == curr->key)
	    return true;
        lr = (item > curr->key);
        prev = curr;
        curr = curr->link[lr];
    }
    return false;
}
 
template <class Type>
inline Type BST<Type>::inOrder(bNode<Type> *ptr) const
{
    bool lr = 1;
    Type temp;
    bNode<Type> *prev = ptr;
 
    ptr = ptr->link[1];
    while (ptr->link[0] != NULL) {
        prev = ptr;
        ptr = ptr->link[lr = 0];
    }
    prev->link[lr] = ptr->link[1];
    temp = ptr->key;
    delete ptr;
    return temp;
}
 
template <class Type>
inline int BST<Type>::subNodes(bNode<Type>* const &ptr) const
{
    if (ptr->link[1] != NULL) {
        if (ptr->link[0] != NULL)
            return 3;
        else
            return 2;
    }
    else if (ptr->link[0] != NULL)
        return 1;
    else
        return 0;
}
 
template <class Type>
int BST<Type>::height(bNode<Type>* const &ptr) const
{
    if (ptr == NULL)
        return 0;
 
    int lt = height(ptr->link[0]), rt = height(ptr->link[1]);
 
    if (lt < rt)
        return 1 + rt;
    return 1 + lt;
}
 
template <class Type>
bNode<Type>* BST<Type>::minmax(bNode<Type> *ptr, const bool &lr) const
{
    while (ptr->link[lr] != NULL)
        ptr = ptr->link[lr];
    return ptr;
}
 
template <class Type>
BST<Type>::BST()
{
    root = NULL;
        count = 0; 
}
 
template <class Type>
BST<Type>::~BST()
{
    clear(root);
}
 
template <class Type>
void BST<Type>::clear()
{
    clear(root);
    root = NULL;
    count = 0; 
}
 
template <class Type>
bool BST<Type>::isEmpty() const
{
    return (root == NULL);
}
 
template <class Type>
bool BST<Type>::insert(const Type &item)
{
    if (root == NULL) {
        root = new bNode<Type>;
        root->key = item;
        root->count++;
        return true;
    }
    bool lr;
    bNode<Type> *curr = root, *prev;
 
    if (search(item, curr, prev, lr))
        return false;
    prev->link[lr] = new bNode<Type>;
    prev->link[lr]->key = item;
    count++;
    return true;
}
 
template <class Type>
bool BST<Type>::remove(const Type &item)
{
    bool lr = 1;
    bNode<Type> *curr = root, *prev;
 
    if (!search(item, curr, prev, lr))
        return false;
    switch (int s = subNodes(curr)) {
    case 0:
    case 1:
    case 2:
        if (curr == root)
            root = curr->link[(s > 1)];
        else
            prev->link[lr] = curr->link[(s > 1)];
        delete curr;
        break;
    case 3:
        curr->key = inOrder(curr);
    }
    count--;
    return true;
}
 
template <class Type>
bool BST<Type>::search(const Type &item, Type *&ptr) const
{
    bool found;
    bNode<Type> *curr = root, *prev;
 
    found = search(item, curr, prev, found);
    ptr = &curr->key;
    return found;
}
 
template <class Type>
Type BST<Type>::min() const
{
    return minmax(root, 0)->key;
}
 
template <class Type>
Type BST<Type>::max() const
{
    return minmax(root, 1)->key;
}
 
template <class Type>
int BST<Type>::size() const
{
    return count;
}
 
template <class Type>
int BST<Type>::height() const
{
    return height(root);
}
