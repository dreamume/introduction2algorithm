#ifndef STRASSEN_H
#define STRASSEN_H

#include <vector>

template <typename T> class Strassen {
public:
    Strassen() = default;
    std::vector<std::vector<T>> compute(const std::vector<std::vector<T>>& a, 
                                        const std::vector<std::vector<T>>& b);
private:
    std::vector<std::vector<T>> computeInternal(const std::vector<std::vector<T>>& a, 
                                                const std::vector<std::vector<T>>& b);
};

#endif
