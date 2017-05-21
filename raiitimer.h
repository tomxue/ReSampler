#ifndef _RAIITIMER_H
#define _RAIITIMER_H 1

#include <iostream>
#include <chrono>

class RaiiTimer {
public:
    RaiiTimer() {
        beginTimer = std::chrono::high_resolution_clock::now();    
    }
    ~RaiiTimer() {
        endTimer = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer - beginTimer).count();
        std::cout << "Time=" << duration << " ms" << std::endl;
    }
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> beginTimer;
    std::chrono::time_point<std::chrono::high_resolution_clock> endTimer;
};

#endif // _RAIITIMER_H