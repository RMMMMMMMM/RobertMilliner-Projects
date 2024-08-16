/*
Robert Michael Milliner
CS 5420
Dr. Rague
Final Project
Due: 04/21/24
Version 1.0
-----------------------------------------------------------------
This program checks if a Digraph has an euler cycle, and if so, finds one

This program takes as input:
    (Optinal): the adjacency list of a Digraph as a '.txt' file

This program produces as output:
    (input given): An euler path for the input graph, if one exists
    (no input): the output for my test cases 
-----------------------------------------------------------------
*/

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <queue>//for BFS
#include "alg_graphs.h"


/*
This is an implementation  of Hierholzer's Algorithm for Digraphs 

Because this algorighm checks every vertex, and traverses every edge exactly once, 
it has time complexity 
    O(v + e)
where v is the number of vertices
and e is the the number of edges

Because I am making a full copy of the Digraph,
it has space complexity
    O(2G)
where G is the entire graph 

This answers my research question! 
Yes, an Euler path in a Digraph can be found in linear time!
*/
std::vector<int> HierholzerDirected(Digraph* g){
    std::vector<int> EulerCycle; // then final euler cylce
    std::vector<int> CurrentPath; //the current path used in the algorithm

    Digraph tempG = Digraph(*g); //temp graph to keep track of what edges have been used 
    int currentVertex = 0; //the vertex we are currently at in the algorithm 
    int nextVertex = 0; // The next vertex to be traversed to in the algorithm

    CurrentPath.push_back(0);//starting at vertex 0

    while(CurrentPath.size() > 0){
        currentVertex = CurrentPath[CurrentPath.size() - 1];

        //the current vertex still has edges in its adj list 
        if(tempG.adj(currentVertex).size() > 0){
            nextVertex = tempG.adj(currentVertex).back();
            tempG.remove_edge(currentVertex, nextVertex);

            CurrentPath.push_back(nextVertex);
        }

        //the current vertex has nothing left in its adj list
        else{
            EulerCycle.push_back(CurrentPath.back());
            CurrentPath.pop_back();
        }
    }

    return EulerCycle;
}

/*
Checks the degree of each vertex in a Digraph

if in-deg == out-deg
    total deg is even for a given vertex

This function has time complexity 
    O(v)
where v is the number of vertices in the graph 
*/
bool checkDegreeDirected(Digraph* g, const int numV){
    for(int i = 0; i < numV; i++){
        if(g->in_degree(i) != g->out_degree(i)){
            return false;
        }
    }
    return true;
}

/*
simple implementation of BFS for directed graphs

I am using this to check for connectivity, so it is returning a bool

BFS will visit each vertex at most once, and each at most once 
thus the time complexity of this alg is:
    O(v + e)
where v is the number of vertices
and e is the number of edges
*/
bool directedBFS(Digraph* g, const int numV, const int startVertex){
    std::queue<int> q;//a queue for BFS
    std::vector<bool> visited(numV, false);//a vector to keep track of what vertices have been visited 

    // Create a copy of the graph to keep track of visited edges
    Digraph tempG = Digraph(*g); 

    visited[startVertex] = true;
    q.push(startVertex);

    while(!q.empty()){
        int currentVertex = q.front();
        q.pop();

        // Add each unvisited vertex that is adjacent to the current vertex to the queue, and mark them as visited, delete the traversed edge
        for(int nextVertex : tempG.adj(currentVertex)){
            if(!visited[nextVertex]){
                visited[nextVertex] = true;
                q.push(nextVertex);

                // Mark the edge as visited
                tempG.remove_edge(currentVertex, nextVertex);
            }
        }
    }

    // Check if all vertices are visited
    for(int i = 0; i < visited.size(); i++){
        if(!visited[i]){
            return false;
        }
    }

    return true;
}

/*
For a directed graph to have a Euler cycle, it must be strongly connected

one way to check for this is to run BFS for every vertex in the graph

because we are running BFS 'v' times, this alg has time complexity
    O(v(v + e)) == O(v^2 + ve)

There are algorithms that could improve this, but the research question for this 
assignment was to investigate the run time of finding Euler Cycles, not strong connectivity
so i went with this approach
*/
bool checkStronglyConnectedDirected(Digraph* g, const int numV){
    for(int i = 0; i < numV; i++){
        if(!directedBFS(g, numV, i)){
            return false;
        }
    }
    return true;
}

/*
For an Euler path to exist in a Digraph:

the in-deg must be equal to the out-deg for each vertex
AND 
the Digraph must be strongly connected 

this function just calls the functions that make the aforementioned checks
and thus has time complexity: 
    O(v^2 + ve + v)
this is just the run times of the two funcitons added together
*/
bool checkExistDirected(Digraph *g, const int numV){
    if(!checkDegreeDirected(g, numV) || !checkStronglyConnectedDirected(g, numV)){
        return false;
    }
    return true;
}

/*
My implementation of Hierholzer's Algorighm stores the euler cycle backwards in a vector

this function just prints a vector backwards and inserts " -> " between each entry

because an euler cycle visits each edge exactly once,
this function has time complexity:
    O(e)
where e is the number of edges in the graph 
*/
void PrintEulerCycle(const std::vector<int>* EulerCycle){    
    for(int i = EulerCycle->size() - 1; i >= 0; i--){
        std::cout << EulerCycle->at(i);
        if(i > 0){
            std::cout << " -> ";
        }
    }
    std::cout << std::endl;
}

/*
takes filename as input

looks for ':'
when it finds ':'
    add edges, FROM current vertex --> TO all V in current list 
*/
Digraph* FileToDiGraph(const std::string filename, int* numV, int* numE){
    //open file
    std::ifstream inputFile(filename);
    //make sure file is actually open
    //if no, terminate program
    if (!inputFile.is_open()) {
        std::cout << "Error opening file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    //temp variables to be used while parsing file
    char tempChar;//tempChar to be used for reading from file
    std::string tempString = "";//tempString to be used for string to int conversion
    bool loopBool = true; //for looping
    int count = 0; //for looping


    //1st line will always be # of V
    while(loopBool){
        inputFile.get(tempChar);
        tempString = tempString + tempChar;//change from char to string
        if(tempChar == '\n'){
            loopBool == false;
            break;
        }
    }
    Digraph* g = new Digraph(std::stoi(tempString));//make digraph with 'tempString' vertices
    *numV = std::stoi(tempString);//set numv
    tempString = "";//reset tempString


    //we need a seperate loop for vertex 0 because of 
    //how i am handling 'count'
    while(loopBool){//vertex 0 
        loopBool = true;
        inputFile.get(tempChar);
        if(tempChar == ':'){
            while(loopBool){
                inputFile.get(tempChar);
                if(tempChar == '\n' || inputFile.eof()){//end of vertex 0 adj list
                    count++;
                    loopBool = false;
                    break;
                }
                else if(tempChar == '0' || tempChar == '1' || tempChar == '2' || 
                        tempChar == '3' || tempChar == '4' || tempChar == '5' || 
                        tempChar == '6' || tempChar == '7' || tempChar == '8' || 
                        tempChar == '9'){//in adj list
                        
                    tempString = tempString + tempChar;//change from char to string

                    if(inputFile.peek() == '0' || inputFile.peek() == '1' || inputFile.peek() == '2' || 
                        inputFile.peek() == '3' || inputFile.peek() == '4' || inputFile.peek() == '5' || 
                        inputFile.peek() == '6' || inputFile.peek() == '7' || inputFile.peek() == '8' || inputFile.peek() == '9'){
                            inputFile.get(tempChar);
                            tempString = tempString + tempChar;//change from char to string
                        }


                    
                    g->add_edge(count, std::stoi(tempString));
                    *numE += 1;//increment edge count
                    tempString = "";//reset tempString
                }
            }
        }
    }
    while(inputFile.get(tempChar)){//all other V 
        loopBool = true;
        if(tempChar == ':'){
            while(loopBool){
                inputFile.get(tempChar);
                if(tempChar == '\n' || inputFile.eof()){//end of current vertex adj list
                    count++;
                    loopBool = false;
                    break;
                }
                else if(tempChar == '0' || tempChar == '1' || tempChar == '2' || 
                        tempChar == '3' || tempChar == '4' || tempChar == '5' || 
                        tempChar == '6' || tempChar == '7' || tempChar == '8' || 
                        tempChar == '9'){//in adj list
                        
                    tempString = tempString + tempChar;//change from char to string

                    if(inputFile.peek() == '0' || inputFile.peek() == '1' || inputFile.peek() == '2' || 
                        inputFile.peek() == '3' || inputFile.peek() == '4' || inputFile.peek() == '5' || 
                        inputFile.peek() == '6' || inputFile.peek() == '7' || inputFile.peek() == '8' || inputFile.peek() == '9'){
                            inputFile.get(tempChar);
                            tempString = tempString + tempChar;//change from char to string
                        }


                    
                    g->add_edge(count, std::stoi(tempString));
                    *numE += 1;//increment edge count
                    tempString = "";//reset tempString
                }
            }
        }
    }
    inputFile.close();
    return g;
}

/*
this function prints the results of a test

the total runtime for any given test is:
    Euler Cycle exists
        O(v^2 + ve + 2v + 2e)
    Euler Cycle does not exist
        O(v^2 + ve + v)

which is interesting! 

*/
void runSingleTest(const std::string filename){
    int numV = 0;//number of vertices in graph
    int numE = 0;//number of edges in graph

    Digraph* g = FileToDiGraph(filename, &numV, &numE);//make digraph
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << filename << ":\n\n";
    if(!checkExistDirected(g, numV)){//no euler cycle
        std::cout << "There does not exist an Euler cycle in this graph\n";
        g->~Digraph();
    }
    else{//euler cycle exists
        std::cout << "There exists an Euler cycle in this graph:\n";
        std::vector<int> EulerCycle = HierholzerDirected(g);
        PrintEulerCycle(&EulerCycle);
        g->~Digraph();
    }
    std::cout << "------------------------------------------------------------------------------\n";
}

/*
this function runs all of my test cases 

the test files must be in the same folder as this program
*/
void runAllTests(){
    runSingleTest("test1.txt");
    runSingleTest("test2.txt");
    runSingleTest("figureDG.txt");
    runSingleTest("test4.txt");
    runSingleTest("test5.txt");
    runSingleTest("test6.txt");
    runSingleTest("test7.txt");
}


int main(int argc, char *argv[]){
    if(argc == 2){//this lets us input any graph 
        std::string filename = argv[1]; //get name of input file from command line arg
        runSingleTest(filename);
    }
    else if(argc > 2){//too many inputs, just run test cases
        std::cout << "-----TOO MANY INPUTS, RUNNING TEST CASES-----\n";
        runAllTests();
    }
    else{//if no graphs are given as input, we run the test cases 
        runAllTests();
    }
    return 0;
}


