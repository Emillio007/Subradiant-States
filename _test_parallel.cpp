#include <iostream>
#include <cmath>
#include <chrono>
#include <ctime>

int main() {

    std::cout.sync_with_stdio(false);

    int cap = pow(2, 20);

    int dim = 2;
    int a[dim][dim];
	a[0][0] = 1;
	a[0][1] = 2;
	a[1][0] = 3;
	a[1][1] = 4;
    int b[dim][dim];
	b[0][0] = 5;
	b[0][1] = 5;
	b[1][0] = 5;
	b[1][1] = 5;
    int c[dim][dim];

    auto start = std::chrono::system_clock::now();
    
    std::cout << "Before loop \n";

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

    std::cout << "After loop \n";

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "Trying to output time \n";

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    std::cout.flush();

    return 0;
}
