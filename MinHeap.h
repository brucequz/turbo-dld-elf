#ifndef MIN_HEAP_H
#define MIN_HEAP_H

#include <vector>

struct DetourObject{
    DetourObject(): originalPathIndex(-1) {};
    double pathMetric;
    double forwardPathMetric;
    int detourStage;
    int startingState;
    int originalPathIndex;         //path that is being detoured from, defaults to -1 to indicate no detours
    
    bool operator<(const DetourObject& obj){
        return pathMetric < obj.pathMetric;
    }
    bool operator>(const DetourObject& obj){
        return pathMetric > obj.pathMetric;
    }
};

class minheap{
public:
    minheap();
    void insert(DetourObject);
    DetourObject pop();
    int size();
private:
    std::vector<DetourObject> detourList;
    void reHeap(int index);
    int parentIndex(int index);
    int rightChildIndex(int index);
    int leftChildIndex(int index);
};


#endif