// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007

#include <iostream>
using namespace std;

int M=3;
int N=3;

int main() {

    int** A = new int*[N];
    for (int i = 0; i < N; i++) {
        A[i] = new int[M];
    }
    A[0][0] = 1;
    A[0][1] = 1;
    A[0][2] = -1;
    A[1][0] = 1;
    A[1][1] = 1;
    A[1][2] = 0;
    A[2][0] = -1;
    A[2][1] = 0;
    A[2][2] = -1;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            cout << "A["<<i<<"]["<<j<<"] = " << A[i][j] << " ";
        }
        cout << endl;
    }

    for (int i = 0; i < N; i++) {
        delete [] A[i];
    }
    delete [] A;

    return 0;
}
