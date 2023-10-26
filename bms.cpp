#include <iostream>
#include <string>
#include <bitset> // bitset to convert to binary
#include <vector>
#include <utility> // pair
#include <cmath> // floor, pow
#include <algorithm> // random shuffle
#include <random> // normal distribution

#define END_SUCCESS 0
#define END_ERROR 1

using namespace std;

/* Function ro print an error message and end with a given exit code */
void printError(string message, int exitCode=END_ERROR) {
    cerr << message << endl;
    exit(exitCode);
}

/* Function to parse arguments and check if they're valid */
string parseArgs(int argc, char **argv) {
    if (argc != 2) {
        printError("Invalid argument! Valid arguments: -d | -e");
    }

    if (string(argv[1]) != string("-e") && string(argv[1]) != string("-d")) {
        printError("Invalid argument! Valid arguments: -d | -e");
    }

    return string(argv[1]);
}

/* Function to check if given char is a 0 or 1 */
bool isOneOrZero(char c) {
    return (c >= 48 && c <= 49);
}

/* Function to check if given character is a digit */
bool isDigit(char c) {
    return (c >= 48 && c <= 57);
}

/* Function to check if given character is a letter */
bool isLetter(char c) {
    return ((c >= 65 && c <= 90) || (c >= 97 && c <= 122));
}

/* Function to remove unwanted characters during decoding */
string makeValidDecoding(string input) {
    string tmp = "";
    for (auto c : input) {
        if (isOneOrZero(c)) {
            tmp += c;
        }
    }
    return tmp;
}

/* Function to remove unwanted characters during encoding */
string makeValidEncoding(string input) {
    string tmp = "";
    for (auto c : input) {
        if (isDigit(c) || isLetter(c)) {
            tmp += c;
        }
    }
    return tmp;
}

/* Function to read whole stdin and return it */
string readInput(string arg) {
    string input;
    string tmp;
    if (arg == "-e") { // if encoding, remove everything but a-z, A-Z, 0-9
        while (getline(cin, tmp)) {
            // test that input characters are valid
            tmp = makeValidEncoding(tmp);
            input += tmp;
        }
    }
    else if (arg == "-d") { // decoding - remove everything but 1 and 0
        while (getline(cin, tmp)) {
            tmp = makeValidDecoding(tmp);
            input += tmp;
        }
    }
    return input;
}

/* Function to create a vector of integers representing the original string in binary */
vector<int> encodeBinary(string ascii) {
    vector<int> numbers;
    for (size_t i = 0; i < ascii.size(); i++) { // convert each ascii char
        auto bit_repre = bitset<8>(ascii[i]).to_string();
        auto first_one = bit_repre.find('1'); // find first number thats not zero
        bit_repre = bit_repre.substr(first_one); // remove zero padding
        for (auto number : bit_repre) {
            numbers.push_back(number-48); // ascii 0 is 48, so 48-48 is 0 a 49-48 is 1, just as I want
        }
    }
    return numbers;
}

/* Function to create a vector of integers from the input */
vector<int> decodeBinary(string ascii) {
    vector<int> numbers;
    for (size_t i = 0; i < ascii.size(); i++) {
        auto number = ascii[i]-48;
        // I dont need to check if its 1 or 0, I already have
        numbers.push_back(number);
    }
    return numbers;
}

void printVector(vector<int> vec, bool pretty=false) {
    if (pretty) {
        for (auto number : vec) {
            cout << number << " ";
        }
    }
    else {
        for (auto number : vec) {
            cout << number;
        }
    }
    cout << endl;
}

void printVector(vector<vector<int>> vec) {
    for(size_t i = 0; i < vec.size(); i++) {
        for (size_t j = 0; j < vec[0].size(); j++) {
            cout << vec[i][j] << " ";
        }   
        cout << endl;
    }
}

void vectorInfo(vector<int> vec) {
    cout << vec.size() << endl;
}

void vectorInfo(vector<vector<int>> vec) {
    cout << vec.size() << " " << vec[0].size() << endl;
}

/* Mozna prepsat na jeden radke inicializce vektoru */
vector<vector<int>> zeroMatrix(size_t rows, size_t cols) {
    vector<vector<int>> block;
    for (size_t i = 0; i < rows; i++) {
        vector<int> row = {};
        for (size_t j = 0; j < cols; j++) {
            row.push_back(0);
        }
        block.push_back(row);
    }
    return block;
}

/* Function to simulate python-esque arr[X:Y] */
void setFirstPartOfVector(vector<vector<int>> &arr, int X, int Y, vector<vector<int>> block) {
    auto ctr = 0;
    for (auto i = X; i < Y; i++) {
        arr[i] = block[ctr++];
    }
}

vector<int> getSlice(vector<vector<int>> arr, int X, int Y, int Z) {
    vector<int> slice;
    for (size_t i = X; (int) i < Y; i++) {
        slice.push_back(arr[i][Z]);
    }
    return slice;
}

/* Fuction to transpose a matrix in vector form */
/* https://stackoverflow.com/a/49445850 */
void transpose(vector<vector<int>> &vec) {

    vector<vector<int>> trans(vec[0].size());

    for (size_t i = 0; i < vec.size(); i++) {
        for (size_t j = 0; j < vec[i].size(); j++) {
            trans[j].push_back(vec[i][j]);
        }
    }
    vec = trans;
}

/** Function to construct the H matrix */
vector<vector<int>> construct_H(size_t codeword, size_t d_c, size_t d_v) {

    auto n_equations = floor((codeword * d_v) / d_c);
    auto block_size = floor(n_equations / d_v);

    size_t rows = floor(n_equations / d_v);
    size_t cols = codeword;

    vector<vector<int>> block = zeroMatrix(rows, cols);
    vector<vector<int>> H = zeroMatrix(n_equations, codeword);

    // Fill the first block with consecutive ones in each block row
    for (auto i = 0; i < block_size; i++) {
        for (auto j = i * d_c; j < (i+1) * d_c; j++) {
            block[i][j] = 1;
        }
    }
    // set the H vector up to the block size to correspond with the set block
    setFirstPartOfVector(H, 0, block_size, block);

    // remaining blocks are permutations of the first block's column
    for (size_t i = 1; i < d_v; i++) {
        transpose(block);
        random_shuffle(block.begin(), block.end());
        transpose(block);
        setFirstPartOfVector(H, i * block_size, (i+1) * block_size, block);
    }
    return H;   
}

size_t largestvalueIndex(vector<int> arr) {
    size_t largest = 0;
    for (size_t i = 1; i < arr.size(); i++) {
        if (arr[i] > arr[i-1]) {
            largest = i;
        }
    }
    return largest;
}

vector<vector<int>> identityVector(size_t n) {
    vector<vector<int>> identity = zeroMatrix(n, n);
    for (size_t i = 0; i < n; i++) {
        identity[i][i] = 1;
    }
    return identity;
}

pair<vector<vector<int>>, vector<vector<int>>> gaussjordan(vector<vector<int>> A, bool toggle) {
    
    auto m = A.size();
    auto n = A[0].size();

    size_t pivot_old = -1;

    vector<vector<int>> P;
    
    if (toggle) {
        P = identityVector(m);
    }

    for (size_t j = 0; j < n; j++) {
        vector<int> filtre_down = getSlice(A, pivot_old+1, (int) m, (int) j);
        size_t pivot = largestvalueIndex(filtre_down) + pivot_old + 1;

        if (A[pivot][j]) {
            pivot_old++;
            if (pivot_old != pivot) {
                auto aux = A[pivot];
                A[pivot] = A[pivot_old];
                A[pivot_old] = aux;

                if (toggle) {
                    aux = P[pivot];
                    P[pivot] = P[pivot_old];
                    P[pivot_old] = aux;
                }
            }
            
            for (size_t i = 0; i < m; i++) {
                if (i != pivot_old && A[i][j]) {
                    if (toggle) {
                        vector<int> tmp = {};
                        for (size_t q = 0; q < P[0].size(); q++) {
                            auto abs_val = abs(P[i][q] - P[pivot_old][q]);
                            tmp.push_back(abs_val);
                        }
                        P[i] = tmp;
                    }
                    vector<int> tmp = {};
                        for (size_t q = 0; q < A[0].size(); q++) {
                            auto abs_val = abs(A[i][q] - A[pivot_old][q]);
                            tmp.push_back(abs_val);
                        }
                    A[i] = tmp;
                }
            }
        }
        if (pivot_old == m-1) {
            break;
        }
    }
    return make_pair(A, P);
}

size_t vecSum(vector<int> arr) {
    size_t sum = 0;
    for (size_t k = 0; k < arr.size(); k++) {
        sum += arr[k];
    }
    return sum;
}

size_t vecSum(vector<vector<int>> arr) {
    size_t sum = 0;
    for (size_t i = 0; i < arr.size(); i++) {
        for (size_t k = 0; k < arr[i].size(); k++) {
            sum += arr[i][k];
        }
    }
    return sum;
}

void moduloMatrix(vector<vector<int>> &arr, int n) {
    for(size_t i = 0; i < arr.size(); i++) {
        for(size_t j = 0; j < arr[i].size(); j++) {
            arr[i][j] = arr[i][j] % n;
        }
    }
}

void moduloMatrix(vector<int> &arr, int n) {
    for(size_t i = 0; i < arr.size(); i++) {
            arr[i] = arr[i] % n;
    }
}

/* Function to multiply matrix by another matrix */
vector<vector<int>> matrixProduct(vector<vector<int>> A, vector<vector<int>> B) {
    if (A[0].size() != B.size()) {
        printError("Matrix multiplication failed");
    }

    vector<vector<int>> multiplied(A.size(), vector<int>(B[0].size(), 0));

    for(size_t i = 0; i < A.size(); i++) {
        for(size_t j = 0; j < B[0].size(); j++) {
            for (size_t k = 0; k < A[0].size(); k++) {
                multiplied[i][j] += (A[i][k] * B[k][j]);
            }
        }
    }
    return multiplied;
}

/* Function to multiply matrix by a vector */
vector<int> matrixProduct(vector<vector<int>> A, vector<int> B) {
    if (A[0].size() != B.size()) {
        printError("Matrix multiplication failed");
    }

    vector<int> multiplied(A.size(), 0);

    for(size_t i = 0; i < A.size(); i++) {
        for (size_t k = 0; k < A[0].size(); k++) {
            multiplied[i] += (A[i][k] * B[k]);
        }
    }
    return multiplied;
}

vector<double> vectorMultiply(vector<double> &A, double val) {
    vector<double> multiplied;
    for (size_t i = 0; i < A.size(); i++) {
        multiplied.push_back(A[i] * val);
    } 
    return multiplied;
}


/* Function to do value^vec for each value in the vector */
vector<int> vectorExponent(vector<int> &vec, int value) {
    vector<int> multiplied = {};
    for (size_t i = 0; i < vec.size(); i++) {
        multiplied.push_back(pow(value, ((int) vec[i])));
    }
    return multiplied;
}

/* Function to add values in 2 vectors */
vector<double> vectorAdd(vector<int> &A, vector<double> &B) {
    vector<double> results;
    if (A.size() != B.size()) {
        printError("Error with encoding during vector addition! They're different size!");
    }
    for (size_t i = 0; i < A.size(); i++) {
        results.push_back(A[i] + B[i]);
    }
    return results;
}

vector<vector<int>> construct_G(vector<vector<int>> H) {

    auto n_code = H[0].size();

    transpose(H);
    auto H_tQ = gaussjordan(H, true);
    transpose(H);

    auto Href_colonnes = H_tQ.first;

    auto tQ = H_tQ.second;

    transpose(Href_colonnes);
    auto Href_diag = gaussjordan(Href_colonnes, false).first;
    transpose(Href_colonnes);

    auto Q = tQ;
    transpose(Q);

    auto n_bits = n_code - vecSum(Href_diag);

    auto Y = zeroMatrix(n_code, n_bits);

    setFirstPartOfVector(Y, n_code-n_bits, (int) Y.size(), identityVector(n_bits));

    auto tG = matrixProduct(Q, Y);
    moduloMatrix(tG, 2);

    return tG;
}

/* Function to create LDPC encoding and decoding matrices */
pair<vector<vector<int>>, vector<vector<int>>> make_ldpc(vector<int> binary_input) {

    size_t codeword = binary_input.size() * 2;
    size_t d_c = binary_input.size();
    size_t d_v = binary_input.size() - 1;

    if (d_v <= 1) {
        printError("Input must be at least 1 character long!");
    }

    vector<vector<int>> H = construct_H(codeword, d_c, d_v);
    vector<vector<int>> G = construct_G(H);

    return make_pair(H, G);
}

/* Function to return a random nubmer from standard normal distribution which has mean 0 and variance 1*/
vector<double> randn(size_t n) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> distr(0.0, 1.0);
    vector<double> ret;
    for (size_t i = 0; i < n; i++) {
        ret.push_back(distr(gen));
    }
    return ret;
}

/** @brief Function to perform the LDPC encoding
 * 
 *  @param tG The coding matrix
 *  @param v Vector of integeres containing the input in binary format
 *  @param snr Signal-To-Noise ratio
 * 
 *  @returns Encoded value
 */
vector<int> encode(vector<vector<int>> tG, vector<int> v, size_t snr) {
    auto d = matrixProduct(tG, v);
    moduloMatrix(d, 2);

    /* 
    //Gaussian noise, just to check it works!
    auto x = vectorExponent(d, -1);
    auto sigma = pow(10, (-snr/20));
    auto rand_ns = randn(x.size());
    auto e = vectorMultiply(rand_ns, sigma);
    auto y = vectorAdd(x, e);
    return y;
    */

    return d;
}

void decode(vector<vector<int>> tG, vector<int> y, size_t snr) {
    printVector(y);
}

int main(int argc, char **argv) {
    string arg = parseArgs(argc, argv);
    string input = readInput(arg);
    vector<int> binary_input;

    if (arg == "-e") {
        binary_input = encodeBinary(input);
    }
    if (arg == "-d") {
        binary_input = decodeBinary(input);
    }

    if (binary_input.size() == 0) {
        printError("Invalid input! Please enter allowed character only!");
    }

    pair<vector<vector<int>>, vector<vector<int>>> H_G = make_ldpc(binary_input);
    size_t snr = 20;

    if (arg == "-e") {
        auto encoded_text = encode(H_G.second, binary_input, snr);
        printVector(encoded_text);
    }
    else if (arg == "-d") {
        /* Pro dekodovani mi staci si ze vstupu zjistit, kolik bitu (jaka byla delka) melo puvodni slovo a podle toho udelam ldpc */
        decode(H_G.first, binary_input, snr);
    }

    return END_SUCCESS;
}
