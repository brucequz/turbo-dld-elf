#include "MinHeap.h"

minheap::minheap() {
  // constructor if necessary
}

int minheap::parentIndex(int index) { return (index - 1) / 2; }
int minheap::leftChildIndex(int index) { return (2 * index + 1); }
int minheap::rightChildIndex(int index) { return (2 * index + 2); }

void minheap::insert(DetourObject detour) {
  detourList.push_back(detour);
  int index = detourList.size() - 1;
  while (index > 0 && detourList[parentIndex(index)] > detourList[index]) {
    std::swap(detourList[parentIndex(index)], detourList[index]);
    index = parentIndex(index);
  }
}

DetourObject minheap::pop() {
  DetourObject detour = detourList[0];
  detourList[0] = detourList[detourList.size() - 1];
  detourList.pop_back();
  reHeap(0);

  return detour;
}

void minheap::reHeap(int index) {
  int leftIndex = leftChildIndex(index);
  int rightIndex = rightChildIndex(index);
  int minDetourIndex = index;
  if (leftIndex < detourList.size() &&
      detourList[leftIndex] < detourList[minDetourIndex])
    minDetourIndex = leftIndex;
  if (rightIndex < detourList.size() &&
      detourList[rightIndex] < detourList[minDetourIndex])
    minDetourIndex = rightIndex;

  if (minDetourIndex != index) {
    std::swap(detourList[index], detourList[minDetourIndex]);
    reHeap(minDetourIndex);
  }
}

int minheap::size() { return detourList.size(); }