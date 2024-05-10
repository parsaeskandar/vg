//
// Created by seeskand on 4/13/24.
// This code is partly based on the code provided here: https://github.com/solangii/b-plus-tree/tree/master
//

#ifndef VG_BPLUSTREE_HPP
#define VG_BPLUSTREE_HPP
bool DEBUG = false;

#include <iostream>

template<typename T>
struct Node {
    bool is_leaf;
    std::size_t degree; // maximum number of children
    std::size_t size; // current number of item
    T *item;
    Node<T> **children;
    Node<T> *parent;
    Node<T> *next;
    Node<T> *prev;

public:
    Node(std::size_t _degree) {// Constructor
        this->is_leaf = false;
        this->degree = _degree;
        this->size = 0;

        T *_item = new T[degree - 1];
        for (int i = 0; i < degree - 1; i++) {
            _item[i] = 0;
        }
        this->item = _item;

        Node<T> **_children = new Node<T> *[degree];
        for (int i = 0; i < degree; i++) {
            _children[i] = nullptr;
        }
        this->children = _children;

        this->parent = nullptr;
        this->next = nullptr; // TODO: check all the instances where we find the next node searching the parent node
        this->prev = nullptr; // TODO: check all the instances where we find the previous node searching the parent node

    }
};

template<typename T>
class BPlusTree {
    Node<T> *root;
    std::size_t degree;

public:
    BPlusTree(std::size_t _degree) {// Constructor
        this->root = nullptr;
        this->degree = _degree;
    }

    ~BPlusTree() { // Destructor
        clear(this->root);
    }

    Node<T> *getroot() {
        return this->root;
    }

    Node<T> *BPlusTreeSearch(Node<T> *node, T key) {
        if (node == nullptr) { // if root is null, return nullptr
            return nullptr;
        } else {
            Node<T> *cursor = node; // cursor finding key

            while (!cursor->is_leaf) { // until cusor pointer arrive leaf
                for (int i = 0; i < cursor->size; i++) { //in this index node, find what we want key
                    if (key < cursor->item[i]) { //find some range, and let find their child also.
                        cursor = cursor->children[i];
                        break;
                    }
                    if (i == (cursor->size) - 1) {
                        cursor = cursor->children[i + 1];
                        break;
                    }
                }
            }

            //search for the key if it exists in leaf node.
            for (int i = 0; i < cursor->size; i++) {
                if (cursor->item[i] == key) {
                    return cursor;
                }
            }

            return nullptr;
        }
    }

    Node<T> *BPlusTreeRangeSearch(Node<T> *node, T key) {
        if (node == nullptr) { // if root is null, return nullptr
            return nullptr;
        } else {
            Node<T> *cursor = node; // cursor finding key

            while (!cursor->is_leaf) { // until cusor pointer arrive leaf
                for (int i = 0; i < cursor->size; i++) { //in this index node, find what we want key
                    if (key < cursor->item[i]) { //find some range, and let find their child also.
                        cursor = cursor->children[i];
                        break;
                    }
                    if (i == (cursor->size) - 1) {
                        cursor = cursor->children[i + 1];
                        break;
                    }
                }
            }
            return cursor;
        }
    }

    int range_search(T start, T end, T *result_data, int arr_length) {
        int index = 0;

        Node<T> *start_node = BPlusTreeRangeSearch(this->root, start);
        Node<T> *cursor = start_node;
        T temp = cursor->item[0];

        while (temp <= end) {
            if (cursor == nullptr) {
                break;
            }
            for (int i = 0; i < cursor->size; i++) {
                temp = cursor->item[i];
                if ((temp >= start) && (temp <= end)) {
                    result_data[index] = temp;
                    index++;
                }
            }
            cursor = cursor->children[cursor->size];
        }
        return index;
    }

    bool search(T data) {  // Return true if the item exists. Return false if it does not.
        return BPlusTreeSearch(this->root, data) != nullptr;
    }

    int find_index(T *arr, T data, int len) {
        int index = 0;
        for (int i = 0; i < len; i++) {
            if (data < arr[i]) {
                index = i;
                break;
            }
            if (i == len - 1) {
                index = len;
                break;
            }
        }
        return index;
    }

    T *item_insert(T *arr, T data, int len) {
        int index = 0;
        for (int i = 0; i < len; i++) {
            if (data < arr[i]) {

                index = i;
                break;
            }
            if (i == len - 1) {
                index = len;
                break;
            }
        }

        for (int i = len; i > index; i--) {
            arr[i] = arr[i - 1];
        }

        arr[index] = data;

        return arr;
    }

    T *run_insert(T *arr, T new_run, size_t &len, int run_length, size_t &inserted_index) {
        size_t index = len;
        for (int i = 0; i < len; i++) {
            if (new_run.start_position < arr[i].start_position) {
                index = i;
                break;
            }
        }
        inserted_index = index;

        if (DEBUG)
            cerr << "++++" << new_run << " " << len << " " << index << "   " << arr[index - 1] << arr[index + 1]
                 << endl;



        // the run should be after a gap run, otherwise it is a redundant run
        if (index != 0 && (arr[index - 1].graph_position.value != 0) && (new_run.graph_position.value != 0)) {
            cerr << "Error: the new run is not inserted correctly. The new run cannot inserted in another run!" << endl;
            return arr;
        }

        if (new_run.graph_position.value == 0) {
//            if (index )
            for (size_t i = len; i > index; i--) {
                arr[i] = arr[i - 1];
            }
            arr[index] = new_run;
            len++;
            return arr;
        }
            // first case is when the new_run doesn't hit previous run or the next run
        else if ((new_run.start_position > arr[index - 1].start_position || index == 0) &&
                 (new_run.start_position + run_length < arr[index].start_position || index == len)) {

            if (DEBUG) cout << "insert run case 1" << endl;
            for (size_t i = len + 1; i > index; i--) {
                arr[i] = arr[i - 2];
            }

            arr[index] = new_run;
            arr[index + 1] = create_gap(new_run.start_position + run_length);
            len += 2;
            inserted_index = index;
        }

            // second case is when we hit the previous run end point and not the next run starting point
        else if (new_run.start_position == arr[index - 1].start_position &&
                 ((new_run.start_position + run_length < arr[index].start_position) || index == len)) {
            if (DEBUG) cout << "insert run case 2" << endl;
            // two cases here, first being the new run graph position not be the same as the previous run graph position
            if (new_run.graph_position.value != arr[index - 2].graph_position.value) {
                if (DEBUG) cout << "insert run case 2.1" << endl;
                // in this case have to remove the gap run and insert the new run in its place and then add a gap for the new run
                arr[index - 1] = new_run;
                inserted_index = index - 1;
                // create a place for the gap run
                for (size_t i = len; i > index; i--) {
                    arr[i] = arr[i - 1];
                }
                arr[index] = create_gap(new_run.start_position + run_length);
                len++;
            }
                // the other case is when the new run graph position is the same as the previous non-gap run graph position
            else if (new_run.graph_position.value == arr[index - 2].graph_position.value) {
                if (DEBUG) cout << "insert run case 2.2" << endl;
                // in this case we just have to change the starting position of the gap in the index-1
                arr[index - 1].start_position += run_length;
                inserted_index = index - 2;
            }

        }
            // third case, when we hit the next run starting point but not the previous one
        else if ((new_run.start_position > arr[index - 1].start_position || index == 0) &&
                 new_run.start_position + run_length == arr[index].start_position) {
            if (DEBUG) cout << "insert run case 3" << endl;
            // two cases here, first being the new run graph_position not be the same as the next run graph_position
            if (new_run.graph_position.value != arr[index].graph_position.value) {
                if (DEBUG) cout << "insert run case 3.1" << endl;
                // in this case we don't need to add a gap after adding the new run
                for (size_t i = len; i > index; i--) {
                    arr[i] = arr[i - 1];
                }
                arr[index] = new_run;
                inserted_index = index;
                len++;
            }
                // the other case is when the new run graph_position is the same as the next non-gap run graph_position
            else if (new_run.graph_position.value == arr[index].graph_position.value) {
                if (DEBUG) cout << "insert run case 3.2" << endl;
                // in this case we just have to change the starting position of the next non-gap run
                arr[index].start_position = new_run.start_position;
                inserted_index = index;
            }
        }
            // forth case, which is when we hit both runs from both sides
        else if (new_run.start_position == arr[index - 1].start_position &&
                 new_run.start_position + run_length == arr[index].start_position) {
            if (DEBUG) cout << "insert run case 4" << endl;
            // two cases here, first being the new run graph_position not be the same as the previous run graph_position
            if ((new_run.graph_position.value != arr[index - 2].graph_position.value) || index == 1) {
                if (DEBUG) cout << "insert run case 4.1" << endl;
                // two cases here too! first being the new run graph_position not be the same as the next run graph_position
                if (new_run.graph_position.value != arr[index].graph_position.value) {
                    if (DEBUG) cout << "insert run case 4.1.1" << endl;

                    // in this case we don't merge anything, but we need to remove the index-1 gap run
                    arr[index - 1] = new_run;
                    inserted_index = index - 1;
                }
                    // the other case is when the new run graph_position is the same as the next non-gap run graph_position
                else if (new_run.graph_position.value == arr[index].graph_position.value) {
                    if (DEBUG) cout << "insert run case 4.1.2" << endl;
                    // in this case we just have to change the starting position of the next non-gap run and remove the index-1 gap run
                    arr[index].start_position = new_run.start_position;
                    for (size_t i = index - 1; i < len - 1; i++) {
                        arr[i] = arr[i + 1];
                    }
                    inserted_index = index - 1;
                    len--;
                }
            }
                // the other case is when the new run graph_position is the same as the previous non-gap run graph_position
            else if (new_run.graph_position.value == arr[index - 2].graph_position.value) {
                if (DEBUG) cout << "insert run case 4.2" << endl;
                // two cases here too! first being the new run graph_position not be the same as the next run graph_position
                if (new_run.graph_position.value != arr[index].graph_position.value) {
                    if (DEBUG) cout << "insert run case 4.2.1" << endl;
                    //1-3 4-0   4-3 6-0   6-4 7-0 //TODO

                    // 1-3 6-4 7-0
                    // in this case we don't merge with the next run, but we need to merge with the previous run
//                    arr[index - 1] = new_run;
                    // just remove the index run
                    for (size_t i = index - 1; i < len - 1; i++) {
                        arr[i] = arr[i + 1];
                    }
                    len--;
                    inserted_index = index - 1;
                }
                    // the other case is when the new run graph_position is the same as the next non-gap run graph_position
                else if (new_run.graph_position.value == arr[index].graph_position.value) {
                    if (DEBUG) cout << "insert run case 4.2.2" << endl;
                    // in this case we have to merge 3 runs together!
                    // in this case we only have to remove the index-1 and index runs
                    for (size_t i = index - 1; i < len - 2; i++) {
                        arr[i] = arr[i + 2];
                    }
                    len -= 2;
                    inserted_index = index - 2;
                }
            }
        } else {
            cerr << "Error: the new run is not inserted correctly" << endl;
        }
        return arr;
    }

    Node<T> **child_insert(Node<T> **child_arr, Node<T> *child, int len, int index) {
        for (int i = len; i > index; i--) {
            child_arr[i] = child_arr[i - 1];
        }
        child_arr[index] = child;
        return child_arr;
    }

    Node<T> *child_item_insert(Node<T> *node, T data, Node<T> *child) {
        int item_index = 0;
        int child_index = 0;
        for (int i = 0; i < node->size; i++) {
            if (data < node->item[i]) {
                item_index = i;
                child_index = i + 1;
                break;
            }
            if (i == node->size - 1) {
                item_index = node->size;
                child_index = node->size + 1;
                break;
            }
        }
        for (int i = node->size; i > item_index; i--) {
            node->item[i] = node->item[i - 1];
        }
        for (int i = node->size + 1; i > child_index; i--) {
            node->children[i] = node->children[i - 1];
        }

        node->item[item_index] = data;
        node->children[child_index] = child;

        return node;
    }

    void InsertPar(Node<T> *par, Node<T> *child, T data) {

        //overflow check
        Node<T> *cursor = par;
        if (cursor->size < this->degree - 1) {//not overflow, just insert in the correct position
            //insert item, child, and reallocate
            cursor = child_item_insert(cursor, data, child);
            cursor->size++;
        } else {//overflow
            //make new node
            auto *Newnode = new Node<T>(this->degree);
            Newnode->parent = cursor->parent;

            //copy item
            T *item_copy = new T[cursor->size + 1];
            for (int i = 0; i < cursor->size; i++) {
                item_copy[i] = cursor->item[i];
            }
            item_copy = item_insert(item_copy, data, cursor->size);

            auto **child_copy = new Node<T> *[cursor->size + 2];
            for (int i = 0; i < cursor->size + 1; i++) {
                child_copy[i] = cursor->children[i];
            }
            child_copy[cursor->size + 1] = nullptr;
            child_copy = child_insert(child_copy, child, cursor->size + 1,
                                      find_index(item_copy, data, cursor->size + 1));

            //split nodes
            cursor->size = (this->degree) / 2;
            if ((this->degree) % 2 == 0) {
                Newnode->size = (this->degree) / 2 - 1;
            } else {
                Newnode->size = (this->degree) / 2;
            }

            for (int i = 0; i < cursor->size; i++) {
                cursor->item[i] = item_copy[i];
                cursor->children[i] = child_copy[i];
            }
            cursor->children[cursor->size] = child_copy[cursor->size];

            for (int i = 0; i < Newnode->size; i++) {
                Newnode->item[i] = item_copy[cursor->size + i + 1];
                Newnode->children[i] = child_copy[cursor->size + i + 1];
                Newnode->children[i]->parent = Newnode;
            }
            Newnode->children[Newnode->size] = child_copy[cursor->size + Newnode->size + 1];
            Newnode->children[Newnode->size]->parent = Newnode;

            T paritem = item_copy[this->degree / 2];

            delete[] item_copy;
            delete[] child_copy;

            //parent check
            if (cursor->parent == nullptr) {//if there are no parent node(root case)
                auto *Newparent = new Node<T>(this->degree);
                cursor->parent = Newparent;
                Newnode->parent = Newparent;

                Newparent->item[0] = paritem;
                Newparent->size++;

                Newparent->children[0] = cursor;
                Newparent->children[1] = Newnode;

                this->root = Newparent;

                //delete Newparent;
            } else {//if there already have parent node
                InsertPar(cursor->parent, Newnode, paritem);
            }
        }
    }


    T create_gap(size_t gap_start_position) {
        T gap = {gap_start_position, 0};
        return gap;
    }

    // handle the merging items between nodes
    void merge_items_next(T *arr, Node<T> *next_node, size_t &len) {
        if (next_node != nullptr) {
            if (next_node->item[0].start_position < arr[len - 1].start_position) {
                cerr << "Error: this run overlaps with the next run!" << endl;
                cerr << "Changing the run length to the maximum it can be before overlapping" << endl;
                len--;
            } else {
                // if we hit the next node
                if (next_node->item[0].start_position == arr[len - 1].start_position) {
                    if (next_node->item[0].graph_position == arr[len - 2].graph_position) {
                        // merging the last item with first item in next node
                        // in this case we have to remove the first item in the next node
                        for (int i = 0; i < next_node->size - 1; i++) {
                            next_node->item[i] = next_node->item[i + 1];
                        }
                        next_node->size--;
                        next_node->children[next_node->size] = next_node->children[next_node->size +
                                                                                   1]; // TODO: check this line
                        next_node->children[next_node->size + 1] = nullptr;
                        len--;
                    } else {
                        // in this case we don't want the gap run, so just decrease the len
                        len--;

                    }

                }
            }
        }
        if (DEBUG) {
            cout << "merging items with the next node output " << endl;
            for (int i = 0; i < len; i++) {
                cout << arr[i] << " ";
            }

        }


    }

    // handle the merging items between nodes
    void merge_items_prev(T *arr, Node<T> *prev_node, size_t &len) {
        if (prev_node != nullptr) {
            if (prev_node->item[prev_node->size - 1].graph_position == arr[0].graph_position) {
                if (DEBUG) cout << "merging items with the previous node" << endl;
                // merging the last item with first item in next node
                // in this case we have to remove the first item in the node)
                for (int i = 0; i < len - 1; i++) {
                    arr[i] = arr[i + 1];
                }
                len--;
            }
        }


    }

    // this function check how can the underflow be handled
    void underflow(Node<T> *cursor) {
        if (DEBUG) {
            cout << "entering underflow" << endl;
            cout << cursor->size << endl;
        }
        bool handled = false;
        if (cursor != this->root) {


            int sib_index = -1; // TODO: maybe change this the pointers cursor->next maybe
            for (int i = 0; i < cursor->parent->size + 1; i++) {
                if (cursor == cursor->parent->children[i]) {
                    sib_index = i;
                }
            }
            int left = sib_index - 1;
            int right = sib_index + 1;

            if (cursor->size < std::ceil(this->degree / 2)) {//underflow case
                if (DEBUG) cout << "underflow" << endl;
                if (left >= 0 && !handled) {// left_sibiling exists
                    Node<T> *leftsibling = cursor->parent->children[left];
                    if (DEBUG) cout << "underflow case 1" << endl;


                    if (leftsibling->size > degree / 2) { //if data number is enough to use this node
                        if (DEBUG) cout << "underflow case 1.1" << endl;
                        T *temp = new T[cursor->size + 1];

                        //copy item
                        for (int i = 0; i < cursor->size; i++) {
                            temp[i] = cursor->item[i];
                        }

                        //insert and rearrange
                        item_insert(temp, leftsibling->item[leftsibling->size - 1], cursor->size);
                        for (int i = 0; i < cursor->size + 1; i++) {
                            cursor->item[i] = temp[i];
                        }
                        cursor->size++;
                        delete[] temp;

                        //pointer edit
                        cursor->children[cursor->size] = cursor->children[cursor->size - 1];
                        cursor->children[cursor->size - 1] = nullptr;

                        //sibling property edit
                        leftsibling->item[leftsibling->size - 1] = 0;
                        leftsibling->size--;
                        leftsibling->children[leftsibling->size] = leftsibling->children[leftsibling->size +
                                                                                         1]; //cursor
                        leftsibling->children[leftsibling->size + 1] = nullptr;

                        //parent property edit
                        cursor->parent->item[left] = cursor->item[0];

                        handled = true;
                    }
                }
                if (right <= cursor->parent->size && !handled) {// right_sibiling exists
                    if (DEBUG) cout << "underflow case 2" << endl;
                    Node<T> *rightsibling = cursor->parent->children[right]; // or the cursor->next

                    if (rightsibling->size > degree / 2) {//if data number is enough to use this node
                        if (DEBUG) cout << "underflow case 2.1" << endl;
                        T *temp = new T[cursor->size + 1];

                        //copy item
                        for (int i = 0; i < cursor->size; i++) {
                            temp[i] = cursor->item[i];
                        }
                        //insert and rearrange
                        item_insert(temp, rightsibling->item[0], cursor->size);
                        for (int i = 0; i < cursor->size + 1; i++) {
                            cursor->item[i] = temp[i];
                        }
                        cursor->size++;
                        delete[] temp;

                        //pointer edit
                        cursor->children[cursor->size] = cursor->children[cursor->size - 1];
                        cursor->children[cursor->size - 1] = nullptr;

                        //sibling property edit
                        for (int i = 0; i < rightsibling->size - 1; i++) {
                            rightsibling->item[i] = rightsibling->item[i + 1];
                        }
                        rightsibling->item[rightsibling->size - 1] = 0;
                        rightsibling->size--;
                        rightsibling->children[rightsibling->size] = rightsibling->children[rightsibling->size +
                                                                                            1]; //cursor
                        rightsibling->children[rightsibling->size + 1] = nullptr;

                        //parent property edit
                        cursor->parent->item[right - 1] = rightsibling->item[0];

                        handled = true;
                    }
                }

                //if sibling is not enough to use their data
                //we have to merge step

                if (left >= 0 && !handled) { // left_sibling exists
                    if (DEBUG) cout << "underflow case 3" << endl;
                    Node<T> *leftsibling = cursor->parent->children[left];


                    //merge two leaf node
                    for (int i = 0; i < cursor->size; i++) {
                        leftsibling->item[leftsibling->size + i] = cursor->item[i];
                    }

                    //edit pointer
                    leftsibling->children[leftsibling->size] = nullptr;
                    leftsibling->size = leftsibling->size + cursor->size;
                    leftsibling->children[leftsibling->size] = cursor->children[cursor->size];
                    leftsibling->next = cursor->next;
                    if (leftsibling->next != nullptr) { leftsibling->next->prev = leftsibling; }

                    //parent property edit
                    Removepar(cursor, left, cursor->parent);

//                    if (cursor != nullptr) {
//                        for (int i = 0; i < cursor->size; i++) {
//                            cursor->item[i] = 0;
//                            cursor->children[i] = nullptr;
//                        }
//                        cout << "here 2 " << endl;
//                        cout << cursor->size << endl;
//                        cursor->children[cursor->size] = nullptr;
//
////                    leftsibling->next = cursor->next;
//
//                        cout << "here 2 " << endl;
//                        delete[] cursor->item;
//                        delete[] cursor->children;
//                        delete cursor;
//
//
//                    }

                    handled = true;

                }
                if (right <= cursor->parent->size && !handled) { // right_sibiling exists
                    if (DEBUG) cout << "underflow case 4" << endl;
                    Node<T> *rightsibling = cursor->parent->children[right];

                    //merge two leaf node
                    for (int i = 0; i < rightsibling->size; i++) {
                        cursor->item[i + cursor->size] = rightsibling->item[i];
                    }
                    //edit pointer
                    cursor->children[cursor->size] = nullptr;
                    cursor->size = rightsibling->size + cursor->size;
                    cursor->children[cursor->size] = rightsibling->children[rightsibling->size];
                    cursor->next = rightsibling->next;
                    if (cursor->next != nullptr) { cursor->next->prev = cursor; }

                    //parent property edit
                    Removepar(rightsibling, right - 1, cursor->parent);

                    for (int i = 0; i < rightsibling->size; i++) {
                        rightsibling->item[i] = 0;
                        rightsibling->children[i] = nullptr;
                    }
                    rightsibling->children[rightsibling->size] = nullptr;

                    delete[] rightsibling->item;
                    delete[] rightsibling->children;
                    delete rightsibling;
                    handled = true;

                }

            } else {
                cerr << "Error: the underflow function called incorrectly! " << endl;
            }

        }
    }


    // this function merge two nodes together, we suppose that all the checks (can get merged) are done before calling this function
    void
    merge_nodes(Node<T> *leftsibling, Node<T> *cursor) { // TODO: what if the parents are not the same for the two nodes

        if (DEBUG) cout << "merging nodes" << endl;
        if (DEBUG) cout << cursor->size << endl;
        int left = -1; // TODO: maybe change this the pointers cursor->next maybe
        for (int i = 0; i < cursor->parent->size + 1; i++) {
            if (cursor == cursor->parent->children[i]) {
                left = i;
                break;
            }
        }



        //merge two leaf node
        for (int i = 0; i < cursor->size; i++) {
            leftsibling->item[leftsibling->size + i] = cursor->item[i];
        }
        //edit pointer
        leftsibling->children[leftsibling->size] = nullptr;
        leftsibling->size = leftsibling->size + cursor->size;
        leftsibling->children[leftsibling->size] = cursor->children[cursor->size];
        leftsibling->next = cursor->next;

        //parent property edit
        Removepar(cursor, left, cursor->parent);
        // print items in cursor node
//        cout << "cursor: ";
//        for (int i = 0; i < cursor->size; i++) {
//            cout << cursor->item[i] << " ";
//        }
//        // print all items in cursor parent node
////        cout << "cursor parent: ";
////        for (int i = 0; i < cursor->parent->size; i++) {
////            cout << cursor->parent->item[i] << " ";
////        }
//        cout << "left: " << left << endl;
//        for (int i = 0; i < cursor->size; i++) {
//            cursor->item[i] = 0;
//            cursor->children[i] = nullptr;
//        }
//        cout << "sure" << endl;
//        cout << cursor->size << endl;
//
//        // print items in cursor node
//        cout << "cursor: ";
//        for (int i = 0; i < cursor->size; i++) {
//            cout << cursor->item[i] << " ";
//        }
//        cursor->children[cursor->size] = nullptr;
//        cout << "sure" << endl;
//        delete[] cursor->item;
//        delete[] cursor->children;
//        delete cursor;


    }

    void insert(T data, size_t run_length) {
        if (this->root == nullptr) { //if the tree is empty
            this->root = new Node<T>(this->degree);
            this->root->is_leaf = true;
            this->root->item[0] = data;
            this->root->size = 1; //
            insert(create_gap(data.start_position + run_length), 0); // TODO: changing the run length to 1?
//            this->root->size++;
        } else { //if the tree has at least one node
            Node<T> *cursor = this->root;
            size_t inserted_index;

            //move to leaf node
            cursor = BPlusTreeRangeSearch(cursor, data); // TODO: I think I don't have to change this function
            if (DEBUG) cout << "adding new item " << data;


            //overflow check
            if (cursor->size < (this->degree - 2)) { // not overflow, just insert in the correct position
                //item insert and rearrange

                cursor->item = run_insert(cursor->item, data, cursor->size, run_length, inserted_index);

                if (inserted_index >= cursor->size - 2 && cursor->next != nullptr) {
                    // merge only if we added some new item at the end of the node
                    merge_items_next(cursor->item, cursor->next, cursor->size);
                }
                if (inserted_index == 0 && cursor->prev != nullptr) {
                    // this means the first item starting position is the same as the new run starting position
                    // this leads to error if the first item of the node is not a gap run
                    // so we know that we already removed the gap run and now we have to check if the prev node last item
                    // graph position is the same as the new run graph position
                    // merge only if we added some new item at the end of the node
                    merge_items_prev(cursor->item, cursor->prev, cursor->size);
                }
                // the cursor or the cursor.next can both underflow here!
//                if (cursor->size < std::ceil(this->degree / 2)) {
//                    cout << "underflow on the cursor node" << endl;
//                    if (cursor->next != nullptr) {
//                        if (cursor->next->size + cursor->size < this->degree) {
//                            cout << "underflow on the next node" << endl;
//                            merge_nodes(cursor, cursor->next);
//                            // print all items in cursor and its parent node
////                            cout << "cursor: ";
////                            for (int i = 0; i < cursor->size; i++) {
////                                cout << cursor->item[i] << " ";
////                            }
////                            cout << '\n' << "cursor parent: ------- ";
////                            for (int i = 0; i < cursor->parent->size; i++) {
////                                cout << cursor->parent->item[i] << " ";
////                            }
//
//                        } else {
//                            underflow(cursor);
//                        }
//                    } else {
//                        underflow(cursor);
//                    }
//                }


                if (cursor->next != nullptr) {
                    if (cursor->next->size < std::ceil(this->degree / 2)) {
                        if (DEBUG) cout << "normal underflow on the next node " << endl;
                        underflow(cursor->next);
                    }
                }
                if (cursor->size < std::ceil(this->degree / 2)) {
                    if (DEBUG) cout << "normal underflow on the cursor node" << endl;
                    underflow(cursor);
                }

                //edit pointer(next node)
                cursor->children[cursor->size] = cursor->children[cursor->size - 1];
                cursor->children[cursor->size - 1] = nullptr;
            } else { //overflow case // TODO: this is not a 100% overflow case, merging the Runs could results in not overflowing

                //copy item
                T *item_copy = new T[cursor->size + 2];
                for (int i = 0; i < cursor->size; i++) {
                    item_copy[i] = cursor->item[i];
                }
                size_t item_copy_size = cursor->size;

                item_copy = run_insert(item_copy, data, item_copy_size, run_length, inserted_index);
                // print all item copy items

                if (inserted_index >= cursor->size - 2 && cursor->next != nullptr) {


                    // merge only if we added some new item at the end of the node
                    if (DEBUG) cout << "merging items" << endl;
                    merge_items_next(item_copy, cursor->next, item_copy_size);
                    // print all items in cursor->next node
                    // check for underflow on the next node
//                    if (cursor->next != nullptr) {
//                        if (cursor->next->size < this->degree / 2){
//                            cout << "underflow on the next node" << endl;
//                            underflow(cursor->next);
//                        }
//                    }
                }
                if (inserted_index == 0 && cursor->prev != nullptr) {
                    if (DEBUG) cout << "merging items prev here" << endl;
                    merge_items_prev(item_copy, cursor->prev, item_copy_size);
                    // print all item_copy items
                }

                // Have to check if this was an overflow case or not

                // did not overflow
                if (item_copy_size < this->degree) {
                    cursor->size = item_copy_size;
                    for (int i = 0; i < item_copy_size; i++) {
                        cursor->item[i] = item_copy[i];
                    }
                } else {

                    // still we might have a non-overflow case
                    // there are two cases which are not overflow cases
                    // 1. the new run and its gap run are inserted at the end of the node
                    // might be the case that the new run or its gap run can get merged into the next node first element

//                    if (inserted_index >= cursor->size - 1) {
//                        if (item_copy_size - this->degree == 0) {
//                            // there is one gap run at the end of the node
//                            Node<T> *next_node = cursor->next;
//                            if (next_node != nullptr) {
//                                // two cases
//
//                                if (next_node->item[0].start_position == item_copy[item_copy_size - 1].start_position) {
//                                    if (next_node->item[0].graph_position == item_copy[item_copy_size -
//                                                                                       2].graph_position) { // merging the last item with first item in next node
//                                        // in this case we have to remove the first item in the next node
//                                        // TODO: The parent pointer is correct, can decide to change the value in the parent or not
//                                        for (int i = 0; i < next_node->size - 1; i++) {
//                                            next_node->item[i] = next_node->item[i + 1];
//                                        }
//                                        next_node->size--;
//                                        next_node->children[next_node->size] = next_node->children[next_node->size +
//                                                                                                   1]; // TODO: check this line
//                                        next_node->children[next_node->size + 1] = nullptr;
//                                        for (int i = 0; i < cursor->size; i++) {
//                                            cursor->item[i] = item_copy[i];
//                                        }
//                                        is_overflow = false;
//                                    } else {
//                                        // in this case we don't want the gap run, so just don't do anything
//                                        is_overflow = false;
//                                        for (int i = 0; i < cursor->size; i++) {
//                                            cursor->item[i] = item_copy[i];
//                                        }
//
//                                    }
//
//                                }
//                            }
//
//                        } else if (item_copy_size - this->degree == 1) { // end of the node special case
//                            // there are two extra runs at the end of the node
//                            Node<T> *next_node = cursor->next;
//                            if (next_node != nullptr) {
//                                // two cases
//                                if (next_node->item[0].start_position == item_copy[item_copy_size - 1].start_position) {
//                                    if (next_node->item[0].graph_position == item_copy[item_copy_size -
//                                                                                       2].graph_position) { // merging the last item with first item in next node
//                                        // in this case we have to remove the first item in the next node
//                                        for (int i = 0; i < next_node->size - 1; i++) {
//                                            next_node->item[i] = next_node->item[i + 1];
//                                        }
//                                        next_node->size--;
//                                        next_node->children[next_node->size] = next_node->children[next_node->size +
//                                                                                                   1]; //TODO: check this line
//                                        next_node->children[next_node->size + 1] = nullptr;
//                                        item_copy_size--;
//                                    } else {
//                                        // in this case we don't want the gap run, so just don't do anything
//                                        item_copy_size--;
//                                    }
//                                }
//                            }
//                        }
//                    }

                    //make new node
                    auto *Newnode = new Node<T>(this->degree);
                    Newnode->is_leaf = true;
                    Newnode->parent = cursor->parent;
                    Newnode->next = cursor->next;

                    //split nodes
                    cursor->size = (item_copy_size) / 2;
                    if ((item_copy_size) % 2 == 0) {
                        Newnode->size = (item_copy_size) / 2;
                    } else {
                        Newnode->size = (item_copy_size) / 2 + 1;
                    }

                    for (int i = 0; i < cursor->size; i++) {
                        cursor->item[i] = item_copy[i];
                    }
                    for (int i = 0; i < Newnode->size; i++) {
                        Newnode->item[i] = item_copy[cursor->size + i];
                    }


                    cursor->children[cursor->size] = Newnode;
                    Newnode->children[Newnode->size] = cursor->children[this->degree - 1];
                    cursor->children[this->degree - 1] = nullptr;
                    if (cursor->next != nullptr) {
                        cursor->next->prev = Newnode;
                    }
                    cursor->next = Newnode;
                    Newnode->prev = cursor;
                    delete[] item_copy;


                    //parent check
                    T paritem = Newnode->item[0];


                    if (cursor->parent == nullptr) {//if there are no parent node(root case)
                        auto *Newparent = new Node<T>(this->degree);
                        cursor->parent = Newparent;
                        Newnode->parent = Newparent;

                        Newparent->item[0] = paritem;
                        Newparent->size++;

                        Newparent->children[0] = cursor;
                        Newparent->children[1] = Newnode;

                        this->root = Newparent;
                    } else {//if there already have parent node
                        if (DEBUG) cout << "inserting to existing parent" << endl;
                        InsertPar(cursor->parent, Newnode, paritem);
                    }

                }
                //TODO: do we need a check for the prev node? I don't think so
                if (cursor->next != nullptr) {
                    if (cursor->next->size < std::ceil(this->degree / 2)) {
                        if (DEBUG) cout << "underflow on the next node" << endl;
                        underflow(cursor->next);
                    }
                }
                if (cursor->size < std::ceil(this->degree / 2)) {
                    if (DEBUG) cout << "underflow on the cursor node" << endl;
                    underflow(cursor);
                }

            }
        }
    }

    void remove(T data) { // Remove an item from the tree.
        //make cursor
        Node<T> *cursor = this->root;

        //move to leaf node
        cursor = BPlusTreeRangeSearch(cursor, data);

        //make sibling index
        int sib_index = -1;
        for (int i = 0; i < cursor->parent->size + 1; i++) {
            if (cursor == cursor->parent->children[i]) {
                sib_index = i;
            }
        }
        int left = sib_index - 1;
        int right = sib_index + 1;


        //find data
        int del_index = -1;
        for (int i = 0; i < cursor->size; i++) {
            if (cursor->item[i] == data) {
                del_index = i;
                break;
            }
        }
        //if data dosen't exist, nothing happen
        if (del_index == -1) {
            return; // there is no match remove value
        }

        //remove data
        for (int i = del_index; i < cursor->size - 1; i++) {
            cursor->item[i] = cursor->item[i + 1];
        }
        cursor->item[cursor->size - 1] = 0;
        cursor->size--;

        //if cursor is root, and there are no more data -> clean!
        if (cursor == this->root && cursor->size == 0) {//root case
            clear(this->root);
            this->root = nullptr;
            return;
        }
        cursor->children[cursor->size] = cursor->children[cursor->size + 1];
        cursor->children[cursor->size + 1] = nullptr;


        //underflow check
        if (cursor == this->root) {
            return;
        }
        if (cursor->size < degree / 2) {//underflow case

            if (left >= 0) {// left_sibiling exists
                Node<T> *leftsibling = cursor->parent->children[left];

                if (leftsibling->size > degree / 2) { //if data number is enough to use this node
                    T *temp = new T[cursor->size + 1];

                    //copy item
                    for (int i = 0; i < cursor->size; i++) {
                        temp[i] = cursor->item[i];
                    }

                    //insert and rearrange
                    item_insert(temp, leftsibling->item[leftsibling->size - 1], cursor->size);
                    for (int i = 0; i < cursor->size + 1; i++) {
                        cursor->item[i] = temp[i];
                    }
                    cursor->size++;
                    delete[] temp;

                    //pointer edit
                    cursor->children[cursor->size] = cursor->children[cursor->size - 1];
                    cursor->children[cursor->size - 1] = nullptr;

                    //sibling property edit
                    leftsibling->item[leftsibling->size - 1] = 0;
                    leftsibling->size--;
                    leftsibling->children[leftsibling->size] = leftsibling->children[leftsibling->size + 1]; //cursor
                    leftsibling->children[leftsibling->size + 1] = nullptr;

                    //parent property edit
                    cursor->parent->item[left] = cursor->item[0];

                    return;
                }
            }
            if (right <= cursor->parent->size) {// right_sibiling exists
                Node<T> *rightsibling = cursor->parent->children[right];

                if (rightsibling->size > degree / 2) {//if data number is enough to use this node
                    T *temp = new T[cursor->size + 1];

                    //copy item
                    for (int i = 0; i < cursor->size; i++) {
                        temp[i] = cursor->item[i];
                    }
                    //insert and rearrange
                    item_insert(temp, rightsibling->item[0], cursor->size);
                    for (int i = 0; i < cursor->size + 1; i++) {
                        cursor->item[i] = temp[i];
                    }
                    cursor->size++;
                    delete[] temp;

                    //pointer edit
                    cursor->children[cursor->size] = cursor->children[cursor->size - 1];
                    cursor->children[cursor->size - 1] = nullptr;

                    //sibling property edit
                    for (int i = 0; i < rightsibling->size - 1; i++) {
                        rightsibling->item[i] = rightsibling->item[i + 1];
                    }
                    rightsibling->item[rightsibling->size - 1] = 0;
                    rightsibling->size--;
                    rightsibling->children[rightsibling->size] = rightsibling->children[rightsibling->size +
                                                                                        1]; //cursor
                    rightsibling->children[rightsibling->size + 1] = nullptr;

                    //parent property edit
                    cursor->parent->item[right - 1] = rightsibling->item[0];

                    return;
                }
            }

            //if sibling is not enough to use their data
            //we have to merge step

            if (left >= 0) { // left_sibling exists
                Node<T> *leftsibling = cursor->parent->children[left];

                //merge two leaf node
                for (int i = 0; i < cursor->size; i++) {
                    leftsibling->item[leftsibling->size + i] = cursor->item[i];
                }
                //edit pointer
                leftsibling->children[leftsibling->size] = nullptr;
                leftsibling->size = leftsibling->size + cursor->size;
                leftsibling->children[leftsibling->size] = cursor->children[cursor->size];

                //parent property edit
                Removepar(cursor, left, cursor->parent);
                for (int i = 0; i < cursor->size; i++) {
                    cursor->item[i] = 0;
                    cursor->children[i] = nullptr;
                }
                cursor->children[cursor->size] = nullptr;

                delete[] cursor->item;
                delete[] cursor->children;
                delete cursor;

                return;

            }
            if (right <= cursor->parent->size) { // right_sibiling exists
                Node<T> *rightsibling = cursor->parent->children[right];

                //merge two leaf node
                for (int i = 0; i < rightsibling->size; i++) {
                    cursor->item[i + cursor->size] = rightsibling->item[i];
                }
                //edit pointer
                cursor->children[cursor->size] = nullptr;
                cursor->size = rightsibling->size + cursor->size;
                cursor->children[cursor->size] = rightsibling->children[rightsibling->size];

                //parent property edit
                Removepar(rightsibling, right - 1, cursor->parent);

                for (int i = 0; i < rightsibling->size; i++) {
                    rightsibling->item[i] = 0;
                    rightsibling->children[i] = nullptr;
                }
                rightsibling->children[rightsibling->size] = nullptr;

                delete[] rightsibling->item;
                delete[] rightsibling->children;
                delete rightsibling;
                return;

            }

        } else {
            return;
        }
    }

    void Removepar(Node<T> *node, int index, Node<T> *par) {
        Node<T> *remover = node;
        Node<T> *cursor = par;
        T target = cursor->item[index];

        if (cursor == this->root && cursor->size == 1) {
            Node<T> *remainingChild = (remover == cursor->children[0]) ? cursor->children[1] : cursor->children[0];

            // Clean up the node to remove
            clear(remover);

            // Set the remaining child as the new root
            this->root = remainingChild;
            this->root->parent = nullptr; // It's important to disconnect the new root from its former parent

            // Clean up the old root node
            delete[] cursor->item;   // Free the array of items
            delete[] cursor->children; // Free the array of child pointers
            delete cursor; // Free the parent node itself
            // print the new root and its parent

            return;

            //if cursor is root, and there are no more data -> child node is to be root!
//        if (cursor == this->root && cursor->size == 1) {//root case
//            if (remover == cursor->children[0]) {
//                cout << "here 1" << endl;
//                delete[] remover->item;
//                delete[] remover->children;
//                delete remover;
//                this->root = cursor->children[1];
//                this->root->parent = nullptr;
//                delete[] cursor->item;
//                delete[] cursor->children;
//                delete cursor;
////                return;
//            }
//            if (remover == cursor->children[1]) {

//                cout << "here" << endl;
//                delete[] remover->item;
//                delete[] remover->children;
//                delete remover;
//
//                this->root = cursor->children[0];
//                this->root->parent = nullptr;
//                // print items in root
//                cout << "root: ";
//                for (int i = 0; i < this->root->size; i++) {
//                    cout << this->root->item[i] << " ";
//                }
//                cout << endl;
//                // print items in cursor node
//                cout << "cursor normal case: ";
//                for (int i = 0; i < cursor->size; i++) {
//                    cout << cursor->item[i] << " ";
//                }
//                //print items in cursor node first children
//                cout << "cursor children: ";
//                for (int i = 0; i < cursor->children[0]->size; i++) {
//                    cout << cursor->children[0]->item[i] << " ";
//                }
//                cout << endl;
//                delete[] cursor->item;
//                delete[] cursor->children;
//                delete cursor;
//                cout << "here" << endl;
//                return;
        } else {


            //remove data
            for (int i = index; i < cursor->size - 1; i++) {
                cursor->item[i] = cursor->item[i + 1];
            }
            cursor->item[cursor->size - 1] = 0;

            //remove pointer
            int rem_index = -1;
            for (int i = 0; i < cursor->size + 1; i++) {
                if (cursor->children[i] == remover) {
                    rem_index = i;
                }
            }
            if (rem_index == -1) {
                return;
            }
            for (int i = rem_index; i < cursor->size; i++) {
                cursor->children[i] = cursor->children[i + 1];
            }
            cursor->children[cursor->size] = nullptr;
            cursor->size--;

            //underflow check
            if (cursor == this->root) {
                return;
            }
            if (cursor->size < degree / 2) {//underflow case

                int sib_index = -1;
                for (int i = 0; i < cursor->parent->size + 1; i++) {
                    if (cursor == cursor->parent->children[i]) {
                        sib_index = i;
                    }
                }
                int left = sib_index - 1;
                int right = sib_index + 1;

                if (left >= 0) {// left_sibiling exists
                    Node<T> *leftsibling = cursor->parent->children[left];

                    if (leftsibling->size > degree / 2) { //if data number is enough to use this node
                        T *temp = new T[cursor->size + 1];

                        //copy item
                        for (int i = 0; i < cursor->size; i++) {
                            temp[i] = cursor->item[i];
                        }

                        //insert and rearrange at cursor
                        item_insert(temp, cursor->parent->item[left], cursor->size);
                        for (int i = 0; i < cursor->size + 1; i++) {
                            cursor->item[i] = temp[i];
                        }
                        cursor->parent->item[left] = leftsibling->item[leftsibling->size - 1];
                        delete[] temp;

                        Node<T> **child_temp = new Node<T> *[cursor->size + 2];
                        //copy child node
                        for (int i = 0; i < cursor->size + 1; i++) {
                            child_temp[i] = cursor->children[i];
                        }
                        //insert and rearrange at child
                        child_insert(child_temp, leftsibling->children[leftsibling->size], cursor->size, 0);

                        for (int i = 0; i < cursor->size + 2; i++) {
                            cursor->children[i] = child_temp[i];
                        }
                        delete[] child_temp;

                        //size edit
                        cursor->size++;
                        leftsibling->size--;
                        return;

                    }
                }

                if (right <= cursor->parent->size) {// right_sibiling exists
                    Node<T> *rightsibling = cursor->parent->children[right];

                    if (rightsibling->size > degree / 2) {//if data number is enough to use this node
                        T *temp = new T[cursor->size + 1];

                        //copy item
                        for (int i = 0; i < cursor->size; i++) {
                            temp[i] = cursor->item[i];
                        }
                        //insert and rearrange at cursor
                        item_insert(temp, cursor->parent->item[sib_index], cursor->size);
                        for (int i = 0; i < cursor->size + 1; i++) {
                            cursor->item[i] = temp[i];
                        }
                        cursor->parent->item[sib_index] = rightsibling->item[0];
                        delete[] temp;

                        //insert and reaarange at child

                        cursor->children[cursor->size + 1] = rightsibling->children[0];
                        for (int i = 0; i < rightsibling->size; i++) {
                            rightsibling->children[i] = rightsibling->children[i + 1];
                        }
                        rightsibling->children[rightsibling->size] = nullptr;

                        cursor->size++;
                        rightsibling->size--;
                        return;

                    }
                }

                //if sibling is not enought to use their data
                //we have to merge step
                if (left >= 0) { // left_sibling exists
                    Node<T> *leftsibling = cursor->parent->children[left];

                    leftsibling->item[leftsibling->size] = cursor->parent->item[left];
                    //merge two leaf node
                    for (int i = 0; i < cursor->size; i++) {
                        leftsibling->item[leftsibling->size + i + 1] = cursor->item[i];
                    }
                    for (int i = 0; i < cursor->size + 1; i++) {
                        leftsibling->children[leftsibling->size + i + 1] = cursor->children[i];
                        cursor->children[i]->parent = leftsibling;
                    }
                    for (int i = 0; i < cursor->size + 1; i++) {
                        cursor->children[i] = nullptr;
                    }
                    leftsibling->size = leftsibling->size + cursor->size + 1;
                    //delete recursion
                    Removepar(cursor, left, cursor->parent);
                    return;

                }
                if (right <= cursor->parent->size) { // right_sibiling exists
                    Node<T> *rightsibling = cursor->parent->children[right];

                    cursor->item[cursor->size] = cursor->parent->item[right - 1];
                    //merge two leaf node
                    for (int i = 0; i < rightsibling->size; i++) {
                        cursor->item[cursor->size + 1 + i] = rightsibling->item[i];
                    }
                    for (int i = 0; i < rightsibling->size + 1; i++) {
                        cursor->children[cursor->size + i + 1] = rightsibling->children[i];
                        rightsibling->children[i]->parent = rightsibling;
                    }
                    for (int i = 0; i < rightsibling->size + 1; i++) {
                        rightsibling->children[i] = nullptr;
                    }
                    //edit pointer
                    rightsibling->size = rightsibling->size + cursor->size + 1;
                    //parent property edit
                    Removepar(rightsibling, right - 1, cursor->parent);
                    return;
                }
            } else {
                return;
            }
        }
    }

    void clear(Node<T> *cursor) {
        if (cursor != nullptr) {
            if (!cursor->is_leaf) {
                for (int i = 0; i <= cursor->size; i++) {
                    clear(cursor->children[i]);
                }
            }
            delete[] cursor->item;
            delete[] cursor->children;
            delete cursor;
        }
    }

    void bpt_print() {
        print(this->root);
    }

    void print(Node<T> *cursor) {
        // You must NOT edit this function.
        if (cursor != NULL) {
            for (int i = 0; i < cursor->size; ++i) {
                std::cout << cursor->item[i] << " ";
            }
            std::cout << "\n";

            if (!cursor->is_leaf) {
                for (int i = 0; i < cursor->size + 1; ++i) {
                    print(cursor->children[i]);
                }
            }
        }
    }

    void bpt_check_items() {
        check_items(this->root);
    }

    void check_items(Node<T> *cursor) {
        if (cursor != NULL) {
            if (cursor->is_leaf) {
                for (int i = 0; i < cursor->size - 1 ; i++) {
                    if (cursor->item[i].graph_position == cursor->item[i + 1].graph_position) {
                        cout << "Error: two adjacent items have the same graph position" << endl;
                        cout << "graph position: " << cursor->item[i].graph_position.value << endl;
                        cout << "start position: " << cursor->item[i].start_position << " " << cursor->item[i+1].start_position << endl;
                    }
                }
            } else {
                for (int i = 0; i < cursor->size + 1; ++i) {
                    check_items(cursor->children[i]);
                }
            }
        }
    }
};

#endif
