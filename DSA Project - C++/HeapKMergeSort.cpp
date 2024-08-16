
#ifndef HEAPKMERGESORT_HPP
#define HEAPKMERGESORT_HPP

template <typename T, unsigned int K>
class MinHeap {
public:
	MinHeap(T* sourceArray, unsigned int sourceArraySize);
	T getSmallest();
private:
	void addSubArray(T* item);
	T* heap[K]; // Make this a stack array for performance. This holds all pointers to existing array values for pointer arithmetic
	unsigned int size{ 0 };
	T* sourceArray{ nullptr };
	unsigned int sourceArraySize{ 0 };
};

template <typename T, unsigned int K>
MinHeap<T, K>::MinHeap(T* sourceArray, unsigned int sourceArraySize) {
	this->sourceArray = sourceArray;
	this->sourceArraySize = sourceArraySize;
	if (sourceArraySize < K) {
        for (unsigned int i = 0; i < sourceArraySize; i++) {
			this->addSubArray(sourceArray + i);
		}	
	} else {
		for (unsigned int i = 0; i < K; i++) {
			this->addSubArray(sourceArray + (i * sourceArraySize) / K);
		}	
	}	
}

template <typename T, unsigned int K>
void MinHeap<T, K>::addSubArray(T* item) {

	unsigned int currIndex = size;
 	heap[size] = item;
	size++;

	// Work upward. Ensure the array starting with the lowest value ends up at the top of the heap tree.
	while (currIndex > 0 && *heap[currIndex] < *heap[(currIndex - 1) / 2]) {
		//swap
		T* temp = heap[currIndex];
		heap[currIndex] = heap[(currIndex - 1) / 2];
		heap[(currIndex - 1) / 2] = temp;
		currIndex = (currIndex - 1) / 2;
	}
}

template <typename T, unsigned int K>
T MinHeap<T, K>::getSmallest() {

	// Get top, replace it with the last item to keep the tree balanced
	T retVal = *heap[0];

	// See if this subarray has more values or not.  
		
	if (((int(heap[0] - sourceArray) + 2) * K - 1) % sourceArraySize >= K) {
		printf("At index %d, sourceArraySize %d, math is %d, bumping one\n", (int(heap[0] - sourceArray)), sourceArraySize, (((int(heap[0] - sourceArray) + 2) * K - 1) % sourceArraySize));
	    if (sourceArraySize == 48) {
			printf("Here");
		}
		// We're still within a subarray, use the next value of the subarray
		heap[0] = heap[0] + 1;
		//return retVal;
	}
	else {
		printf("At index %d sourceArraySize %d, math is %d, --- subarray done\n", (int(heap[0] - sourceArray)), sourceArraySize, (((int(heap[0] - sourceArray) + 2) * K - 1) % sourceArraySize));
	
		// This sub array is done.  Put the last leaf of the heap at root
		heap[0] = heap[size - 1];
		size--;
	}

	// Work downward
	unsigned int parent = 0;
	do {
		unsigned int left = 2 * parent + 1;
		unsigned int right = 2 * parent + 2;
		if (left >= size) {
			return retVal;
		}
		if (right < size) {
			if (*heap[left] < *heap[parent]) {
				if (*heap[right] < *heap[left]) {
					T* temp = heap[parent];
					heap[parent] = heap[right];
					heap[right] = temp;
					parent = right;
				}
				else {
					T* temp = heap[parent];
					heap[parent] = heap[left];
					heap[left] = temp;
					parent = left;
				}
			}
			else if (*heap[right] < *heap[parent]) {
				T* temp = heap[parent];
				heap[parent] = heap[right];
				heap[right] = temp;
				parent = right;
			}
			else {
				return retVal;
			}
		}
		else {
			if (*heap[left] < *heap[parent]) {
				T* temp = heap[parent];
				heap[parent] = heap[left];
				heap[left] = temp;
			}
			return retVal;
		}
	} while (true);

}

template <typename T, unsigned int K>
class HeapKMerge : public BaseSort<T> {
public:
	HeapKMerge(const unsigned int capacity) : BaseSort<T>("Heap - K-Way Merge", capacity) {}
	void runSort();
private:
	void runSort(unsigned int firstIndex, unsigned int lastIndex);

};

template <typename T, unsigned int  K>
void HeapKMerge<T, K>::runSort() {
	runSort(0, this->capacity);
}

template <typename T, unsigned int K>
void HeapKMerge<T, K>::runSort(unsigned int firstIndex, unsigned int lastIndex) {

	if (lastIndex - firstIndex < 2) {
		return;
	}

	unsigned int numValues = lastIndex - firstIndex;
	for (unsigned int i = 0; i < K; i++) {
		runSort(firstIndex + numValues * i / K, firstIndex + numValues * (i + 1) / K);
	}

	T* arrayCopy = new T[numValues];
	for (unsigned int i = 0; i < numValues; i++) {
		arrayCopy[i] = this->arr[firstIndex + i];
	}

	// Build heap.  First put all values in the heap.
	MinHeap<T, K> minHeap(arrayCopy, numValues);

	unsigned int itemsSorted = 0;

	while (itemsSorted < numValues) {

		this->arr[firstIndex + itemsSorted] = minHeap.getSmallest();

		itemsSorted++;
	}
	delete[] arrayCopy;

}

#endif