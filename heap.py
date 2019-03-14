import heapq

class Heap:
    def __init__(self, lt_comp=lambda x, y: x < y):
        self.heap = []
        self.lt_comp = lt_comp

    def push(self, element):
        class Element:
            def __init__(self, value, lt_comp):
                self.lt_comp = lt_comp
                self.value = value

            def __lt__(self, other):
                return self.lt_comp(self.value, other.value)

        heapq.heappush(self.heap, Element(element, self.lt_comp))

    def pop(self):
        return heapq.heappop(self.heap).value

    def peek(self):
        return self.heap[0].value

    def peek_bottom(self, functor):
        largest = None
        for ele in self.heap:
            if functor(ele.value):
                if largest is None or largest < ele:
                    largest = ele
        return largest.value

    def empty(self):
        return len(self.heap) == 0

    def delete(self, functor):
        delete_idxs = []
        for idx, ele in enumerate(self.heap):
            if functor(ele.value):
                delete_idxs.append(idx)
        # reverse order so that we delete from the last to the first element
        delete_idxs.reverse()
        for idx in delete_idxs:
            del self.heap[idx]
        heapq.heapify(self.heap)

    def __len__(self):
        return len(self.heap)

    def __iter__(self):
        for ele in self.heap:
            yield ele.value