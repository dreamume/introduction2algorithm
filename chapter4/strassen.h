#ifndef STRASSEN_H
#define STRASSEN_H

#include <vector>

class Strassen {
public:
    Strassen() = default;
    std::vector<std::vector<int>> compute(const std::vector<std::vector<int>>& a, 
                                          const std::vector<std::vector<int>>& b);
private:
    std::vector<std::vector<int>> computeInternal(const std::vector<std::vector<int>>& a, 
                                                  const std::vector<std::vector<int>>& b);
};

#endif
