#include <iostream>
#include <cmath>
#include <chrono>
#include <ctime>

int main() {

    int cap = pow(2, 20);

    int dim = 2;
    int a[dim][dim] = {{1, 2}, {3, 4}};
    int b[dim][dim] = {{5, 5}, {5, 5}};
    int c[dim][dim];

    auto start = std::chrono::system_clock::now();

    for(int l = 0; l < cap; l++){    
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                c[i][j] = 0;
                for (int k = 0; k < dim; k++) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
    }

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}