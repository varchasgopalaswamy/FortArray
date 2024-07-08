#include <iostream>
#include "fortarray.h"
#include <chrono>

int main() {
    const int nx = 2000, ny = 2000;
    const int n = nx*ny;
    auto a = Array::Array<float, 2>(nx,ny);
    auto avec = a.Buffer();
    std::cout << a.Size() << std::endl;
    std::cout << a.State() << std::endl;

    double end;
    auto start = std::chrono::steady_clock::now();
    for (int k = 0 ; k < 100; k++) {
        for (int i = 0; i < a.Shape(0); i++) {
            for (int j = 0; j < a.Shape(1); j++) {
                //const size_t idx = i*a.Shape(1) + j;
                //avec[idx] = static_cast<float>(i + j);
                a[i,j] = static_cast<float>(i + j);
            }
        }
    }
    end = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "Time: " << end << std::endl;

    start = std::chrono::steady_clock::now();
    auto b = a + 1.0;
    end = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "Time: " << end << std::endl;
    std::cout << a[0,0] << ", " << b[0,0] << std::endl;

    start = std::chrono::steady_clock::now();
    auto c = b;
    c += 2.0;
    end = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "Time: " << end << std::endl;
    std::cout << a[0,0] << ", " << b[0,0] << ", " << c[0,0] << std::endl;

    start = std::chrono::steady_clock::now();
    auto d = a + b;
    end = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "Time: " << end << std::endl;
    std::cout << a[0,0] << ", " << b[0,0] << ", " << d[0,0] << std::endl;

    start = std::chrono::steady_clock::now();
    auto e = a * b;
    end = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "Time: " << end << std::endl;
    std::cout << a[0,0] << ", " << b[0,0] << ", " << d[0,0] << std::endl;

    // float b[n];
    // for (int i = 0; i < a.Shape(0); i++) {
    //     for (int j = 0; j < a.Shape(1); j++) {
    //         const size_t idx = i*a.Shape(1) + j;
    //         avec[idx] = static_cast<float>(i + j);
    //         b[i*nx + j] = a[i,j];
    //     }
    // }

    // for (int i = 0; i < n; i++) {
    //     std::cout << avec[i] << " ";
    // }
    auto size = {0};

    auto f = d.View({0});
    return 0;
}
