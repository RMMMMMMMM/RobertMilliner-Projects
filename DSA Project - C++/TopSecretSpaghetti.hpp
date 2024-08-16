#include <iostream>
#include <string>
#include <utility>

using std::cout;
using std::endl;
using std::string;
using std::pair;


/********************************************
******************* Nodes *******************
*********************************************/
template<typename T>
struct SLLNode{//Node with one pointer
public:
    T data{};
    SLLNode<T>* next{nullptr};
};

template<typename T>
struct DLLNode{//Node with two pointers 
public:
    T data{};
    DLLNode<T>* next{nullptr};
    DLLNode<T>* prev{nullptr};
};

template<typename K, typename V>
struct KVNode{//Node with Key and Value
public:
    K key{};
    V value{};
    KVNode<K, V>* next{nullptr};
};

/********************************************
************ Singly-Linked-List *************
*********************************************/
template<typename T>
class SinglyLinkedList{
public:
    void insert(const T& data);//inserts from the front
    void insertBack(const T& data);
    void remove();//delete from front
    void deleteBack();
    void update(const unsigned int index, const T& data);
    T get(const unsigned int index);
    string getName();

private:
    SLLNode<T>* head{nullptr};
    SLLNode<T>* tail{nullptr};
};


template<typename T>
void SinglyLinkedList<T>::insert(const T& data){
    SLLNode<T>* temp = new SLLNode<T>();
    temp->data = data;

    if(!this->head){//empty list
        this->head = temp;
        this->tail = temp;
    }
    else{//one or more nodes
        temp->next = this->head;
        this->head = temp;
    }
}

template<typename T>
void SinglyLinkedList<T>::insertBack(const T& data){
    SLLNode<T>* temp = new SLLNode<T>();
    temp->data = data;

    if(!this->head){//empty list
        this->head = temp;
        this->tail = temp;
    }
    else{//one or more nodes
        SLLNode<T>* curr = this->head;
        while(curr->next){
            curr = curr->next;
        }
        curr->next = temp;
        this->tail = temp;
    }
}

template<typename T>
void SinglyLinkedList<T>::remove(){
    if(!this->head){//empty list
        cout << "The list is already empty" << endl;
    }
    else if(this->head == this->tail){//one node
        delete this->head;
        this->head = nullptr;
        this->tail = nullptr;
        cout << "The list is now empty" << endl;
    }
    else{//more than one node
        SLLNode<T>* curr = this->head;
        this->head = this->head->next;
        delete curr;
    }
}

template<typename T>
void SinglyLinkedList<T>::deleteBack(){
    if(!this->head){//empty list
        cout << "The list is already empty" << endl;
    }
    else if(this->head == this->tail){//one node
        delete this->head;
        this->head = nullptr;
        this->tail = nullptr;
        cout << "The list is now empty" << endl;
    }
    else{//more than one node
        SLLNode<T>* curr = this->head;
        while(curr->next->next){
            curr = curr->next;
        }
        delete curr->next;
        curr->next = nullptr;
        this->tail = curr;
    }
}

template<typename T>
void SinglyLinkedList<T>::update(const unsigned int index, const T& data){//indexing starts at 0
    if(!this->head){//empty list
        cout << "The list is empty..." << endl;
    }
    else{//one or more nodes
        SLLNode<T>* curr = this->head;
        for(unsigned int i = 0; i < index; i++){
            if(curr->next){
                curr = curr->next;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                return;
            }
        }
        curr->data = data;
    }


    // //this code does NOT check boundary cases 
    // SLLNode<T>* curr = this->head;
    // for(unsigned int i = 0; i < index; i++){
    //     if(curr->next){
    //         curr = curr->next;
    //     }
    // }
    // curr->data = data;


}

template<typename T>
T SinglyLinkedList<T>::get(const unsigned int index){
    if(!this->head){
        cout << "The list is empty.." << endl;
        throw std::out_of_range("Empty List");
    }
    else{
        SLLNode<T>* curr = this->head;
        for(unsigned int i = 0; i < index; i++){
            if(curr->next){
                curr = curr->next;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                throw std::out_of_range("Out of bounds");
            }
        }
        return curr->data;
    }



    // //this code does NOT check boundary cases 
    // SLLNode<T>* curr = this->head;
    // for(unsigned int i = 0; i < index; i++){
    //     if(curr->next){
    //         curr = curr->next;
    //     }
    // }
    // return curr->data;
    
}

template<typename T>
string SinglyLinkedList<T>::getName(){
    return "Singly-Linked-List";
}

/********************************************
************ Doubly-Linked-List *************
*********************************************/
template<typename T>
class DoublyLinkedList{
public:
    void insert(const T& data);//inserts from the front 
    void insertBack(const T& data);
    void remove();//delete from front
    void deleteBack();
    void updateFront(const unsigned int index, const T& data);
    void updateBack(const unsigned int index, const T& data);
    void update(const unsigned int index, const T& data);//uses front or back depending on index
    T getFront(const unsigned int index);
    T getBack(const unsigned int index);
    T get(const unsigned int index);//uses front or back depending on index
    string getName();

private:
    DLLNode<T>* head{nullptr};
    DLLNode<T>* tail{nullptr};
    unsigned int count = 0;//only used for update/get methods
};


template<typename T>
void DoublyLinkedList<T>::insert(const T& data){
    DLLNode<T>* temp = new DLLNode<T>();
    temp->data = data;

    if(!this->head){//empty list
        this->head = temp;
        this->tail = temp;
    }
    else{//one or more nodes
        temp->next = this->head;
        this->head->prev = temp;
        this->head = temp;
    }
    this->count++;
}

template<typename T>
void DoublyLinkedList<T>::insertBack(const T& data){
    DLLNode<T>* temp = new DLLNode<T>();
    temp->data = data;

    if(!this->head){//empty list
        this->head = temp;
        this->tail = temp;
    }
    else{//one or more nodes
        temp->prev = this->tail;
        this->tail->next = temp;
        this->tail = temp;
    }
    this->count++;
}

template<typename T>
void DoublyLinkedList<T>::remove(){
    if(!this->head){//empty list
        cout << "The list is already empty" << endl;
    }
    else if(this->head == this->tail){//one node
        delete this->head;
        this->head = nullptr;
        this->tail = nullptr;
        cout << "The list is now empty" << endl;
        this->count--;
    }
    else{//more than one node
        DLLNode<T>* curr = this->head;
        this->head = this->head->next;
        this->head->prev = nullptr;
        delete curr;
        this->count--;
    }
}

template<typename T>
void DoublyLinkedList<T>::deleteBack(){
    if(!this->head){//empty list
        cout << "The list is already empty" << endl;
    }
    else if(this->head == this->tail){//one node
        delete this->head;
        this->head = nullptr;
        this->tail = nullptr;
        cout << "The list is now empty" << endl;
        this->count--;
    }
    else{//more than one node
        DLLNode<T>* curr = this->tail;
        this->tail = this->tail->prev;
        this->tail->next = nullptr;
        delete curr;
        this->count--;
    }
}

template<typename T>
void DoublyLinkedList<T>::updateFront(const unsigned int index, const T& data){
    if(!this->head){//empty list
        cout << "The list is empty..." << endl;
    }
    else{//one or more nodes
        DLLNode<T>* curr = this->head;
        for(unsigned int i = 0; i < index; i++){
            if(curr->next){
                curr = curr->next;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                return;
            }
        }
        curr->data = data;
    }
}

template<typename T>
void DoublyLinkedList<T>::updateBack(const unsigned int index, const T& data){
    if(!this->head){//empty list
        cout << "The list is empty..." << endl;
    }
    else{//one or more nodes
        DLLNode<T>* curr = this->tail;
        for(unsigned int i = 0; i < index; i++){
            if(curr->prev){
                curr = curr->prev;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                return;
            }
        }
        curr->data = data;
    }
}

template<typename T>
void DoublyLinkedList<T>::update(const unsigned int index, const T& data){
    if(!this->head){//empty list
        cout << "The list is empty..." << endl;
    }
    else if(index < (this->count / 2)){
        //updateFront(index, data);
    
        DLLNode<T>* curr = this->head;
        for(unsigned int i = 0; i < index; i++){
            if(curr->next){
                curr = curr->next;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                return;
            }
        }
        curr->data = data;
    }
    else{
        //updateBack(index, data);

        DLLNode<T>* curr = this->tail;
        for(unsigned int i = 0; i < index; i++){
            if(curr->prev){
                curr = curr->prev;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                return;
            }
        }
        curr->data = data;
    }
}

template<typename T>
T DoublyLinkedList<T>::getFront(const unsigned int index){
    if(!this->head){
        cout << "The list is empty.." << endl;
        throw std::out_of_range("Empty List");
    }
    else{
        DLLNode<T>* curr = this->head;
        for(unsigned int i = 0; i < index; i++){
            if(curr->next){
                curr = curr->next;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                throw std::out_of_range("Out of bounds");
            }
        }
        return curr->data;
    }
}

template<typename T>
T DoublyLinkedList<T>::getBack(const unsigned int index){
    if(!this->head){
        cout << "The list is empty.." << endl;
        throw std::out_of_range("Empty List");
    }
    else{
        DLLNode<T>* curr = this->tail;
        for(unsigned int i = 0; i < index; i++){
            if(curr->prev){
                curr = curr->prev;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                throw std::out_of_range("Out of bounds");
            }
        }
        return curr->data;
    }
}

template<typename T>
T DoublyLinkedList<T>::get(const unsigned int index){
    if(index < (this->count / 2)){
        return getFront(index);
    }
    else{
        return getBack(index);
    }
}

template<typename T>
string DoublyLinkedList<T>::getName(){
    return "Doubly-Linked-List";
}

/********************************************
*************** Array-Queue *****************
*********************************************/
template<typename T>
class ArrayQueue{
public:
    ArrayQueue(const unsigned int capacity);
    void insert(const T& data);//enqueue
    void remove();//dequeue
    void update(const unsigned int index, const T& data);
    T get(const unsigned int index); //does this make sense for queues?
    string getName();

private:
    unsigned int head = 0;
    unsigned int tail = 0;
    unsigned int capacity = 0;
    unsigned int size = 0;
    T* arr{nullptr};
};


    template<typename T>
    ArrayQueue<T>::ArrayQueue(const unsigned int capacity){
        this->capacity = capacity;
        this->arr = new T[capacity];
    }

    template<typename T>
    void ArrayQueue<T>::insert(const T& data){
        if(this->size >= this->capacity){//full queue
            cout << "The queue is already full..." << endl;
            return;
        }
        if(this->tail == this->capacity){
            this->tail = 0;
        }
        this->arr[tail] = data;
        this->tail++;
        this->size++;
    }

    template<typename T>
    void ArrayQueue<T>::remove(){
        if(!size){
            cout << "The queue is already empty..." << endl;
            return;
        }
        this->head++;
        if(this->head == this->capacity){
            this->head == 0;
        }
        this->size--;
    }

    template<typename T>
    void ArrayQueue<T>::update(const unsigned int index, const T& data){
        // //This code lets you insert outside of the queue
        // //
        // //  EX: |head|1|2|3|4|5|tail| | |
        // //  index:    0 1 2 3 4      5 6
        // //  ArrayQueue(6,6)
        // //      |head|1|2|3|4|5|tail| |6|
        // // 
        // //index 6 was updated but it is outside of the queue
        //
        // if(index > this->capacity){
        //     cout << "Out of bounds" << endl;
        //     return;
        // }
        // else{
        //     this->arr[index] = data;
        // }


        //this code fixes the aforementioned issue
        if(index > this->capacity){
            cout << "The queue is empty..." << endl;
            return;
        }
        else if(this->head < this-> tail){
            if(index < this->head || index > this->tail){
                cout << "Out of bounds" << endl;
                return;
            }
        }
        else if(this->head > this-> tail){
            if(index < this->tail || index > this->head){
                cout << "Out of bounds" << endl;
                return;
            }
        }
        else{//head == tail
            if(index != this->head){
                cout << "Out of bounds" << endl;
                return;
            }
        }
        this->arr[index] = data;
    }

    template<typename T>
    T ArrayQueue<T>::get(const unsigned int index){
        if(index > this->capacity){
            cout << "The queue is empty..." << endl;
            throw std::out_of_range("Empty Queue");
        }
        else if(this->head < this-> tail){
            if(index < this->head || index > this->tail){
                cout << "Out of bounds" << endl;
                throw std::out_of_range("Out of bounds");
            }
        }
        else if(this->head > this-> tail){
            if(index < this->tail || index > this->head){
                cout << "Out of bounds" << endl;
                throw std::out_of_range("Out of bounds");
            }
        }
        else if(this->head == this->tail){
            if(index != this->head){
                cout << "Out of bounds" << endl;
                throw std::out_of_range("Out of bounds");
            }
        }
        return this->arr[index];
    }

    template<typename T>
    string ArrayQueue<T>::getName(){
        return "Array-Queue";
    }

/********************************************
**************** List-Queue *****************
*********************************************/
template<typename T>
class ListQueue{
public:
    //these functions are identical to SLL should I just call those methods instead?
    void insert(const T& data);//enqueue
    void remove();//dequeue
    void update(const unsigned int index, const T& data);
    T get(const unsigned int index);
    string getName();

private:
    SLLNode<T>* head{nullptr};
    SLLNode<T>* tail{nullptr};
};


    template<typename T>
    void ListQueue<T>::insert(const T& data){
        SLLNode<T>* temp = new SLLNode<T>();
        temp->data = data;

        if(!this->head){//empty queue
            this->head = temp;
            this->tail = temp;
        }
        else{//one or more nodes
            SLLNode<T>* curr = this->head;
            while(curr->next){
                curr = curr->next;
            }
            curr->next = temp;
            this->tail = temp;
    }
    }

    template<typename T>
    void ListQueue<T>::remove(){
        if(!this->head){//empty queue
            cout << "The queue is already empty" << endl;
        }
        else if(this->head == this->tail){//one node
            delete this->head;
            this->head = nullptr;
            this->tail = nullptr;
            cout << "The queue is now empty" << endl;
        }
        else{//more than one node
            SLLNode<T>* curr = this->head;
            this->head = this->head->next;
            delete curr;
        }
    }

    template<typename T>
    void ListQueue<T>::update(const unsigned int index, const T& data){
        if(!this->head){//empty queue
            cout << "The queue is empty..." << endl;
        }
        else{//one or more nodes
            SLLNode<T>* curr = this->head;
            for(unsigned int i = 0; i < index; i++){
                if(curr->next){
                    curr = curr->next;
                }
                else{
                    cout << "OUT OF BOUNDS" << endl;
                    return;
                }
            }
        curr->data = data;
        }
    }

template<typename T>
T ListQueue<T>::get(const unsigned int index){
    if(!this->head){
        cout << "The queue is empty.." << endl;
        throw std::out_of_range("Empty List");
    }
    else{
        SLLNode<T>* curr = this->head;
        for(unsigned int i = 0; i < index; i++){
            if(curr->next){
                curr = curr->next;
            }
            else{
                cout << "OUT OF BOUNDS" << endl;
                throw std::out_of_range("Out of bounds");
            }
        }
    return curr->data;
    }
}

template<typename T>
string ListQueue<T>::getName(){
    return "List-Queue";
}

/********************************************
****************** Stack ********************
*********************************************/
template<typename T>
class Stack{
public:
    Stack(const unsigned int capacity);
    void insert(const T& data);
    void remove();
    void update(const unsigned int index, const T& data);
    T get(const unsigned int index);
    string getName();

private:
    unsigned int capacity = 0;
    unsigned int index = 0;
    T* arr{nullptr};
};


    template<typename T>
    Stack<T>::Stack(const unsigned int capacity){
        this->capacity = capacity;
        this->arr = new T[capacity];
    }

    template<typename T>
    void Stack<T>::insert(const T& data){
        if(this->index >= this->capacity){
            cout << "The stack is already full..." << endl;
            return;
        }
        this->arr[this->index] = data;
        this->index++;
    }

    template<typename T>
    void Stack<T>::remove(){
        if(!this->index){
            cout << "The stack is already empty..." << endl;
            return;
        }
        this->index--;
    }

    template<typename T>
    void Stack<T>::update(const unsigned int index, const T& data){
        if(index >= this->index){
            cout << "Out of bounds" << endl;
            return;
        }
        this->arr[index] = data;
    }

    template<typename T>
    T Stack<T>::get(const unsigned int index){
        if(index >= this->index){
            cout << "Out of bounds" << endl;
            throw std::out_of_range("Out of bounds");
        }
        return this->arr[index];
    }
    
template<typename T>
string Stack<T>::getName(){
    return "Stack";
}

/********************************************
****************** Array ********************
*********************************************/
template<typename T>
class Array{
public:
    Array(const unsigned int capacity);
    void insert(const unsigned int index, T& value);//insert and update seem like the same things with arrays 
    void remove(const unsigned int index);
    T get(const unsigned int index);
    string getName();

private:
    unsigned int capacity = 0;
    T* arr{nullptr};
};


template<typename T>
Array<T>::Array(const unsigned int capacity){
    this->capacity = capacity;
    this->arr = new T[capacity];
}

template<typename T>
void Array<T>::insert(const unsigned int index, T& value){//indexing starts at 0
    if(index > this->capacity - 1){
        cout << "OUT OF BOUNDS" << endl;
    }
    else{
        this->arr[index] = value;
    }
}

template<typename T>
void Array<T>::remove(const unsigned int index){
    if(index > this->capacity - 1){
        cout << "OUT OF BOUNDS" << endl;
    }
    else{
        this->arr[index] = nullptr;//idk if this is correct 
    }
}

template<typename T>
T Array<T>::get(const unsigned int index){
    return this->arr[index];
}

/********************************************
****************** KVList *******************
*********************************************/
template<typename K, typename V>
class KVList{//Singly-Linked-List that uses Key-Value-Nodes(KVNodes), used for open hashing 
public:
    void insert(const K& key, const V& value);//inserts from the front
    void remove(const K& key);
    void update(const K& key, const V& value);
    V get(const K& key);
    string getName();//Dont think this will be used but everything else has it ¯\_(ツ)_/¯

private:
    KVNode<K, V>* head{nullptr};
    KVNode<K, V>* tail{nullptr};
};


template<typename K, typename V>
void KVList<K, V>::insert(const K& key, const V& value){
    KVNode<K, V>* temp = new KVNode<K, V>();
    temp->key = key;
    temp->value = value;

    if(!this->head){//empty list
        this->head = temp;
        this->tail = temp;
    }
    else{//one or more nodes
        temp->next = this->head;
        this->head = temp;
    }
}

template<typename K, typename V>
void KVList<K, V>::remove(const K& key){
    if(!this->head){//empty list
        cout << "The list is already empty" << endl;
        return;
    }


    //which 'one node' case is faster?
    //assuming its the 1st one cause we only need to compare keys if head == tail
    //made same assumption in update() and get()
    else if(this->head == this->tail){//one node
        if(this->head->key == key){
            delete this->head;
            this->head = nullptr;
            this->tail = nullptr;
            cout << "The list is now empty" << endl;
            return;
        }
    }
    // else if(this->head == this->tail && this->head->key == key){//one node
    //     delete this->head;
    //     this->head = nullptr;
    //     this->tail = nullptr;
    //     cout << "The list is now empty" << endl;
    //     return;
    // }




    else{//more than one node
        KVNode<K, V>* curr = this->head->next;
        KVNode<K, V>* currPrev = this->head;
        while(curr->next){
            if(curr->key == key){
                currPrev->next = curr->next;
                delete curr;
                return;
            }
            currPrev = curr;
            curr = curr->next;
        }
    }
    cout << "Key Not Found in list" << endl;
}

template<typename K, typename V>
void KVList<K, V>::update(const K& key, const V& value){//indexing starts at 0
    if(!this->head){//empty list
        cout << "The list is empty" << endl;
    }
    else if(this->head == this->tail){//one node
        if(this->head->key == key){
        this->head->key = key;//pretty sure this line isnt needed 
        this->head->value = value;
        return;
        }
    }
    else{//one or more nodes
        KVNode<K, V>* curr = this->headt;
        while(curr->next){
            if(curr->key == key){
                curr->key = key;//pretty sure this line isnt needed 
                curr->value = value;
                return;
            }
            curr = curr->next;
        }
    }
    cout << "Key Not Found in list" << endl;
}

template<typename K, typename V>
V KVList<K, V>::get(const K& key){
    if(!this->head){
        cout << "The list is empty.." << endl;
        throw std::out_of_range("Empty List");
    }
    else if(this->head == this->tail){//one node
        if(this->head->key == key){
            return this->head->value;
        }
        else{
            cout << "Key Not Found in list" << endl;
            throw std::out_of_range("Key Not Found in list");
        }
    }
    else{//one or more nodes 
        KVNode<K, V>* curr = this->headt;
        while(curr->next){
            if(curr->key == key){
                return curr->value;
            }
            curr = curr->next;
        }
        cout << "Key Not Found in list" << endl;
        throw std::out_of_range("Key Not Found in list");
    }
}

template<typename K, typename V>
string KVList<K, V>::getName(){
    return "KVList";
}

/********************************************
************** Open-Hash-Table **************
*********************************************/
template<typename K, typename V>
class OpenHashTable{
public:
    OpenHashTable(const unsigned int numBuckets);
    unsigned int hash(const K& key);
    void insert(const K& key, const V& value);
    void remove(const K& key);
    void update(const K& key, const V& value);
    V get(const K& key);
    string getName();
private:
    unsigned int numBuckets = 10;
    //KVList<K, V> table[];
    DoublyLinkedList<std::pair<K, V>>* table{nullptr};
};

    template<typename K, typename V>
    OpenHashTable<K, V>::OpenHashTable(const unsigned int numBuckets){
        this->numBuckets = numBuckets;
        this->table = new DoublyLinkedList<pair<K, V>>[this->numBuckets]; 
        // for(int i = 0; i < this->numBuckets; i++){
        //     this->table[i] = new DoublyLinkedList<pair<K, V>>;
        // }
    }

    template<typename K, typename V>
    unsigned int OpenHashTable<K, V>::hash(const K& key){
        return key % this->numBuckets;
    }

    template<typename K, typename V>
    void OpenHashTable<K, V>::insert(const K& key, const V& value){
        this->table[hash(key)].insert(pair<K, V>(key, value));
    }

    template<typename K, typename V>
    void OpenHashTable<K, V>::remove(const K& key){
        this->table[hash(key)].remove(key);
    }

    template<typename K, typename V>
    void OpenHashTable<K, V>::update(const K& key, const V& value){
        this->table[hash(key)].update(pair<K, V>(key, value));
    }

    template<typename K, typename V>
    V OpenHashTable<K, V>::get(const K& key){
        return this->table[hash(key)].get(key);
    }

    template<typename K, typename V>
    string OpenHashTable<K, V>::getName(){
        return "Open-Hash-Table";
    }

// /********************************************
// ************* Closed-Hash-Table *************
// *********************************************/
// template <typename K, typename V>
// class ClosedHashTable{
// public:
//     //should this have refactor method? or just included in insert?
//     ClosedHashTable(const unsigned int numBuckets);
//     void insert(const K& key, const V& value);
//     void remove(const K& key);
//     void update(const K& key, const V& value);
//     V get(const K& key);
//     string getName();

// private:
//     unsigned int numBuckets = 10;
//     unsigned int size = 0;
//     pair<K, V>* table{nullptr};
//     int* statusArray{nullptr};
// }; 


//     template<typename K, typename V>
//     ClosedHashTable<K, V>::ClosedHashTable(const unsigned int numBuckets){
//         this->numBuckets = numBuckets;
//         this->statusArray = new int[numBuckets];
//         for(int i = 0; i < numBuckets; i++){
//             this->statusArray[i] = 0;
//         }
//         this->table = new pair<K, V>[numBuckets];
//     }

//     template<typename K, typename V>
//     void ClosedHashTable<K, V>::insert(const K& key, const V& value){
//         // does this only work for ints?
//         unsigned int hash = key % this->numBuckets;

//         while(hash < this->numBuckets){
//             if(this->statusArray[hash] < 1){
//                 this->table.first = key;
//                 this->table.second = value;
//                 this->statusArray[hash] = 1;
//             }
//             count++;
//             hash = (hash + 1) % this->numBuckets;
//         }
//     }

//     template<typename K, typename V>
//     void ClosedHashTable<K, V>::remove(const K& key){

//     }




//     template<typename K, typename V>
//     string ClosedHashTable<K, V>::getName(){
//         return "Closed-Hash-Table";
//     }