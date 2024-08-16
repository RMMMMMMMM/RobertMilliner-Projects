#include "TopSecretSpaghetti.hpp"
#include <iostream>
#include <random>
#include <chrono>
#include <string>

using std::cin;
using std::cout;
using std::endl;
using std::string;



/*
    1. Make hash-tables with all CRUD operations
    
    2. Make Graphs with operations discussed in discord

    3. Make 
        -Copy constructor 
        -Copy assignment
        -Destructor
        -Move constructor
        -Move assignment
       For each Data Structure  

    g++ TopSecretSpaghetti.cpp -o TopSecretSpaghetti.x -03
    ./TopSecretSpaghetti.x
*/


template <typename CONTAINER_TYPE>
float insertTest(CONTAINER_TYPE& container, const unsigned int size) {
    cout << "Beginning " << container.getName() << " Insertion..." << endl; 
    auto t1 = std::chrono::high_resolution_clock::now(); //start timer 
    for(unsigned int i = 0; i < size; i++){
        container.insert(i);
    }
    auto t2 = std::chrono::high_resolution_clock::now(); //end timer
    cout << "Complete..." << endl;
    std::chrono::duration<double, std::milli> fp_ms = t2-t1;
    return fp_ms.count();
}

template <typename CONTAINER_TYPE>
float removeTest(CONTAINER_TYPE& container, const unsigned int size) {
    cout << "Beginning " << container.getName() << " Deletion..." << endl; 
    auto t1 = std::chrono::high_resolution_clock::now(); //start timer 
    for(unsigned int i = 0; i < size; i++){
        container.remove();
    }
    auto t2 = std::chrono::high_resolution_clock::now(); //end timer
    cout << "Complete..." << endl;
    std::chrono::duration<double, std::milli> fp_ms = t2-t1;
    return fp_ms.count();
}

template <typename CONTAINER_TYPE>
float updateTest(CONTAINER_TYPE& container, const unsigned int size) {
    cout << "Beginning " << container.getName() << " Update..." << endl; 
    auto t1 = std::chrono::high_resolution_clock::now(); //start timer 
    for(unsigned int i = 0; i < size; i++){
        container.update(i, i * 2);//(i * 2) so the value actually changes
    }
    auto t2 = std::chrono::high_resolution_clock::now(); //end timer
    cout << "Complete..." << endl;
    std::chrono::duration<double, std::milli> fp_ms = t2-t1;
    return fp_ms.count();
}

template <typename CONTAINER_TYPE>
float getTest(CONTAINER_TYPE& container, const unsigned int size) {
    cout << "Beginning " << container.getName() << " Read..." << endl; 
    auto t1 = std::chrono::high_resolution_clock::now(); //start timer 
    for(unsigned int i = 0; i < size; i++){
        container.get(i);
    }
    auto t2 = std::chrono::high_resolution_clock::now(); //end timer
    cout << "Complete..." << endl;
    std::chrono::duration<double, std::milli> fp_ms = t2-t1;
    return fp_ms.count();
}

int main(){
    
    
    //initialization of random number stuff 
    std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distrib(0, 10000);
    int randomNumber = 0;
    //initialization of timer stuff
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    float wallTime{ 0.0 };




    unsigned int size = 10'000;
    cout << "----------------------------------------------------" << endl;
    cout << size << " Elements" << endl;
    cout << "----------------------------------------------------" << endl;
    
    // //Make each container type
    SinglyLinkedList<int> SLLContainer; 
    DoublyLinkedList<int> DLLContainer;
    ListQueue<int> LQContainer;
    ArrayQueue<int> AQContainer(size);
    Stack<int> StackContainer(size);
    

    // //insert all
    wallTime = insertTest(SLLContainer, size);
    cout << wallTime << endl;
    wallTime = insertTest(DLLContainer, size);
    cout << wallTime << endl;
    wallTime = insertTest(LQContainer, size);
    cout << wallTime << endl;
    wallTime = insertTest(AQContainer, size);
    cout << wallTime << endl;
    wallTime = insertTest(StackContainer, size);
    cout << wallTime << endl;

    cout << "----------------------------------------------------" << endl;
    cout << "----------------------------------------------------" << endl;

    //update all
    wallTime = updateTest(SLLContainer, size);
    cout << wallTime << endl;
    wallTime = updateTest(DLLContainer, size);
    cout << wallTime << endl;
    wallTime = updateTest(LQContainer, size);
    cout << wallTime << endl;
    wallTime = updateTest(AQContainer, size);
    cout << wallTime << endl;
    wallTime = updateTest(StackContainer, size);
    cout << wallTime << endl;

    cout << "----------------------------------------------------" << endl;
    cout << "----------------------------------------------------" << endl;

    //read all
    //needs try/catch for error handling 
    wallTime = getTest(SLLContainer, size);
    cout << wallTime << endl;
    wallTime = getTest(DLLContainer, size);
    cout << wallTime << endl;
    wallTime = getTest(LQContainer, size);
    cout << wallTime << endl;
    wallTime = getTest(AQContainer, size);
    cout << wallTime << endl;
    wallTime = getTest(StackContainer, size);
    cout << wallTime << endl;

    cout << "----------------------------------------------------" << endl;
    cout << "----------------------------------------------------" << endl;

    //delete all
    wallTime = removeTest(SLLContainer, size);
    cout << wallTime << endl;
    wallTime = removeTest(DLLContainer, size);
    cout << wallTime << endl;
    wallTime = removeTest(LQContainer, size);
    cout << wallTime << endl;
    wallTime = removeTest(AQContainer, size);
    cout << wallTime << endl;
    wallTime = removeTest(StackContainer, size);
    cout << wallTime << endl;

    cout << "----------------------------------------------------" << endl;
    cout << "----------------------------------------------------" << endl;








    /*
        OPEN HASH TEST
            cant use existing functions cause need to pass in key and value
    */
    // OpenHashTable<int, int> OpenHashContainer(size);


    // cout << "Beginning " << OpenHashContainer.getName() << " Insertion..." << endl; 
    // t1 = std::chrono::high_resolution_clock::now(); //start timer 
    // for(unsigned int i = 0; i < size; i++){
    //     OpenHashContainer.insert(i, i + 1);
    // }
    // t2 = std::chrono::high_resolution_clock::now(); //end timer
    // cout << "Complete..." << endl;
    // cout << fp_ms.count() << endl;

    // cout << "Beginning " << OpenHashContainer.getName() << " Update..." << endl; 
    // t1 = std::chrono::high_resolution_clock::now(); //start timer 
    // for(unsigned int i = 0; i < size; i++){
    //     OpenHashContainer.update(i, i + 1);
    // }
    // t2 = std::chrono::high_resolution_clock::now(); //end timer
    // cout << "Complete..." << endl;
    // cout << fp_ms.count() << endl;




    return 0;
}