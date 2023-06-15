//The declaration and definition of a tempalte MUST be in the SAME file
#ifndef SIMPLEQ_H
#define SIMPLEQ_H
#include<memory.h>
#include<limits>

namespace ns3 {
    //template declaration
    template<typename eType = double>
    class cirQue {
    public:
        cirQue(int n);
        void PutElem(eType elem);
        eType GetMean();
        eType GetMin();
        void Reset();

    private:
        eType sum;
        const int totalCount;
        int curCount;
        int curPos;
        eType* elemArray;
    };

    //template definition
    template<typename eType>
    cirQue<eType>::cirQue(int n)
        : totalCount(n), curCount(0), curPos(-1), sum(static_cast<eType>(0)) {
        elemArray = new eType[n]();
    }

    template<typename eType>
    void cirQue<eType>::PutElem(eType elem) {
        curPos = (++curPos == totalCount) ? 0 : curPos;
        sum -= elemArray[curPos];
        elemArray[curPos] = elem;
        sum += elem;
        if (curCount < totalCount)
            curCount++;
    }

    template<typename eType>
    eType cirQue<eType>::GetMean() {
        return sum / curCount;
    }

    template<typename eType>
    void cirQue<eType>::Reset() {
        curCount = 0;
        curPos = -1;
        sum = static_cast<eType>(0);
        memset(elemArray, 0, totalCount * sizeof(eType));
    }

    template<typename eType>
    eType cirQue<eType>::GetMin() {
        eType minval = std::numeric_limits<eType>::max();
        for (int i = 0;i < curCount;i++)
            if (elemArray[i] < minval)
                minval = elemArray[i];
        return minval;
    }
}
#endif
