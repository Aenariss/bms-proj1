/**
 * BMS Project
 * Author: Vojtech Fiala <xfiala61> 
 * Note: Very heavily inspired by https://github.com/hichamjanati/pyldpc
*/


#include <iostream>
#include <fstream>
#include <string>
#include <bitset> // bitset to convert to binary
#include <vector>
#include <utility> // pair
#include <cmath> // floor, pow
#include <algorithm> // random shuffle
#include <random> // normal distribution

#define END_SUCCESS 0
#define END_ERROR 1
#define MATRIX_FILE "matica.csv"
#define ITERATION_LIMIT 500

using namespace std;

/* Function ro print an error message and end with a given exit code */
void printError(string message, int exitCode=END_ERROR) {
    cerr << message << endl;
    exit(exitCode);
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


void moduloMatrix(vector<vector<int>> &arr, int n) {
    for(size_t i = 0; i < arr.size(); i++) {
        for(size_t j = 0; j < arr[i].size(); j++) {
            arr[i][j] = arr[i][j] % n;
        }
    }
}

void moduloMatrix(vector<int> &arr, int n) {
    for(size_t i = 0; i < arr.size(); i++) {
            arr[i] = (n + (arr[i] % n)) % n;
    }
}

/**
 * Function to calculate binary product of 2 matrices
*/
vector<int> binaryProduct(vector<vector<int>> H, vector<int> x) {
    auto prod = matrixProduct(H, x);
    moduloMatrix(prod, 2);
    return prod;
}

/**
 * Function to calculate binary product of 2 matrices
*/
int binaryProduct(vector<int> H, vector<int> x) {
    int prod = 0;
    if (H.size() != x.size()) {
        printError("Invalid binary product calculation!");
    }
    for (size_t i = 0; i < H.size(); i++) {
        prod += H[i] * x[i];
    }
    return (2 + (prod % 2)) % 2;
}

/* Function to parse arguments and check if they're valid */
pair<string, string> parseArgs(int argc, char **argv) {
    if (argc != 2 && argc != 4) {
        printError("Invalid argument! Valid arguments: -d | -e | -m");
    }

    string e_d_flag = "";
    string matrix_file = "";

    for (int i = 0; i < argc; i++) {
        auto curr_arg = string(argv[i]);
        if (curr_arg == "-e" || curr_arg == "-d") {
            e_d_flag = curr_arg;
        }
        // Unnecessary matrix was used
        if (curr_arg == "-m") {
            // There has to be atleast 1 more argument which I assume to be the matrix file
            if (i+1 < argc) {
                matrix_file = argv[i+1];
                i++;
            }
            else {
                printError("When using the -m parameter, you need to also give a file with the matrix!");
            }
        }
    }

    if (matrix_file == "" && e_d_flag == "-d") {
        printError("When using the -d parameter, you need to also give the matrix!");
    }

    if (e_d_flag == "") {
        printError("Missing -e or -d!");
    }

    return make_pair(e_d_flag, matrix_file);
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
    string input = "";
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
            break;
        }
    }
    return input;
}

void printVector(vector<int> vec, bool pretty=true) {
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

/* Function to create a vector of integers representing the original string in binary */
vector<int> encodeBinary(string ascii) {
    vector<int> numbers;
    for (size_t i = 0; i < ascii.size(); i++) { // convert each ascii char
        auto bit_repre = bitset<8>(ascii[i]).to_string();
        //auto first_one = bit_repre.find('1'); // find first number thats not zero
        //bit_repre = bit_repre.substr(first_one); // remove zero padding
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

vector<int> getSlice(vector<int> arr, int X, int Y) {
    vector<int> slice;
    for (size_t i = X; (int) i < Y; i++) {
        slice.push_back(arr[i]);
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
    int largestVal = arr[0];
    for (size_t i = 1; i < arr.size(); i++) {
        if (arr[i] > largestVal) {
            largest = i;
            largestVal = arr[i];
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

vector<double> vectorMultiply(vector<double> &A, double val) {
    vector<double> multiplied;
    for (size_t i = 0; i < A.size(); i++) {
        multiplied.push_back(A[i] * val);
    } 
    return multiplied;
}

vector<int> vectorMultiply(vector<int> &A, double val) {
    vector<int> multiplied;
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
vector<int> encode(vector<vector<int>> tG, vector<int> v) {
    auto d = binaryProduct(tG, v);
    return d;
}

/**
 * Function to simulate numpy where
*/
vector<vector<int>> where(vector<vector<int>> H) {
    vector<int> row;
    vector<int> col;
    for (size_t i = 0; i < H.size(); i++) {
        for (size_t k = 0; k < H[i].size(); k++) {
            if(H[i][k]) {
                col.push_back(k);
                row.push_back(i);
            }
        }
    }
    vector<vector<int>> rowCols = {row, col};
    return rowCols;
}

/***
 * Function to count number of occurences of a given number in a given array
*/
int count(vector<int> arr, int x) {
    int ctr = 0;
    for (auto number : arr) {
        if (number == x) {
            ctr++;
        }
    }
    return ctr;
}

/**
 * Function to simulate the binCount numpy function
*/
vector<int> binCount(vector<int> arr) {
    vector<int> bins;

    int max = *max_element(arr.begin(), arr.end());

    for (int i = 0; i < max+1; i++) {
        auto n_of_occurences = count(arr, i);
        bins.push_back(n_of_occurences);
    }
    return bins;
}

/**
 * Function to find the bits and nodes of the parity matrix H
*/
vector<vector<int>> bitsNodes(vector<vector<int>> H) {

    auto validIndexes = where(H);
    auto bitsIndices = validIndexes[0];
    auto bits = validIndexes[1];

    transpose(H);

    validIndexes = where(H);
    auto nodesIndices = validIndexes[0];
    auto nodes = validIndexes[1];

    auto bitsBins = binCount(bitsIndices);
    auto nodeBins = binCount(nodesIndices);

    vector<vector<int>> bitsNodes = {bitsBins, bits, nodeBins, nodes};
    return bitsNodes;
}

/**
 * Function to calculate the probabilities
*/
pair<pair<vector<vector<int>>, vector<vector<int>>>, vector<double>> log_belief_propagation(vector<int> bitsHist, vector<int> bits, 
                                                                          vector<int> nodesHist, vector<int> nodes, 
                                                                          vector<vector<int>> Lq, vector<vector<int>> Lr, vector<int> Lc, int n_iter) {

    auto m = Lr.size();
    auto n = Lr[0].size();
    auto n_messages = 1;
    auto bits_counter = 0;
    auto nodes_counter = 0;

    // Horizontal
    for (size_t i = 0; i < m; i++) {
        auto ff = bitsHist[i];
        auto ni = getSlice(bits, bits_counter, bits_counter + ff);
        bits_counter += ff;
        for (auto j : ni) {
            auto nij = ni;
            vector<int> X = {1}; // number of n_messages
            if (n_iter == 0) {
                for (size_t k = 0; k < nij.size(); k++) {
                    if (nij[k] != j) {
                        for (size_t kk = 0; kk < X.size(); kk++) {
                            X[kk] *= tanh(0.5 * Lc[nij[kk]]);
                        }
                    }
                }
            }
            else {
                for (size_t k = 0; k < nij.size(); k++) {
                    if (nij[k] != j) {
                        for (size_t kk = 0; kk < X.size(); kk++) {
                            X[kk] *= tanh(0.5 * Lq[i][nij[kk]]);
                        }
                    }
                }
            }
            auto num = X;
            for (size_t kk = 0; kk < num.size(); kk++) {
                num[kk] += 1;
            }
            auto denom = X;
            for (size_t kk = 0; kk < denom.size(); kk++) {
                denom[kk] = 1 - denom[kk];
            }
            for (int kk = 0; kk < n_messages; kk++) {
                if (num[kk] == 0) {
                    Lr[i][j] = -1;
                }
                else if (denom[kk] == 0) {
                    Lr[i][j] = 1;
                }
                else {
                    Lr[i][j] = log(num[kk] / denom[kk]);
                }
            }
        }
    }

    // Vertical
    for (size_t j = 0; j < n; j++) {
        auto ff = nodesHist[j];
        auto mj = getSlice(nodes, nodes_counter, nodes_counter + ff);
        nodes_counter += ff;
        for (auto i : mj) {
            auto mji = mj;
            Lq[i][j] = Lc[j];
            for (size_t k = 0; k < mji.size(); k++) {
                if (mji[k] != i) {
                    Lq[i][j] += Lr[mji[k]][j];
                }
            }
        }
    }

    // Posterior
    vector<double> L_Post;
    for (size_t i = 0; i < n; i++) {
        L_Post.push_back(0);
    }
    nodes_counter = 0;

    for (size_t j = 0; j < n; j++) {
        auto ff = nodesHist[j];
        auto mj = getSlice(nodes, nodes_counter, nodes_counter + ff);
        nodes_counter += ff;
        double value = 0;
        for (auto number : mj) {
            value += Lr[number][j];
        }
        L_Post[j] = Lc[j] + value;
    }

    return make_pair(make_pair(Lq, Lr), L_Post);
}

int checkMatrixZero(vector<vector<int>> H, vector<int> x) {
    auto prod = binaryProduct(H, x);
    for (auto a : prod) {
        if (a != 0) {
            return 0;
        }
    }
    return 1;
}

pair<vector<vector<int>>, vector<int>> gaussElimination(vector<vector<int>> A, vector<int> b) {
    auto n = A.size();
    auto k = A[0].size();

    auto range_limit = n < k ? n : k;

    for (size_t j = 0; j < range_limit; j++) {
        vector<int> listedpivots;
        for (size_t i = j; i < n; i++) {
            if (A[i][j]) {
                listedpivots.push_back(i);
            }
        }
        size_t pivot = 99999999;
        if (listedpivots.size()) {
            pivot = *min_element(listedpivots.begin(), listedpivots.end());
        }
        else {
            continue;
        }

        if (pivot != j) {
            auto aux = A[j];
            A[j] = A[pivot];
            A[pivot] = aux;

            auto baux = b[j];
            b[j] = b[pivot];
            b[pivot] = baux;
        }

        for (auto i = j+1; i < n; i++) {
            if (A[i][j]) {
                for (size_t k = 0; k < A[0].size(); k++) {
                    A[i][k] = abs(A[i][k]-A[j][k]);
                }
                b[i] = abs(b[i]-b[j]);
            }
        }
    }
    return make_pair(A,b);
}

void invertBits(vector<int>&message) {
    for (size_t i = 0; i < message.size(); i++) {
        message[i] = message[i] == 0 ? 1 : 0;
    }
}

void binaryToString(vector<int> &message) {

    string tmp = "";

    for (size_t i = 0; i < message.size(); i++) {
        if (i % 8 == 0 && i > 0) {
            auto ascii_val = stoi(tmp, nullptr, 2);
            cout << (char) ascii_val;
            tmp = to_string(message[i]);
        }
        else {
            tmp += to_string(message[i]);
        }
    }
    auto ascii_val = stoi(tmp, nullptr, 2);
    cout << (char) ascii_val; // Finish the last one
    cout << endl;
}

/**
 * Function to get the original message from the decoded sequence and print it
*/
void originalMessage(vector<vector<int>> H, vector<int> x) {
    vector<vector<int>> G = construct_G(H);

    auto k = G[0].size();
    auto gaussResult = gaussElimination(G, x);
    auto rtG = gaussResult.first;
    auto rx = gaussResult.second;

    vector<int> message;
    for (size_t i = 0; i < k; i++) {
        message.push_back(0);
    }

    message[k-1] = rx[k-1];
    for (int i = k-1-1; i >= 0; i--) {
        message[i] = rx[i];

        vector<int> rtgSlice;
        vector<int> messageSlice;
        for (size_t z = i+1; z < k; z++) {
            rtgSlice.push_back(rtG[i][z]);
            messageSlice.push_back(message[z]);
        }
        message[i] -= binaryProduct(rtgSlice, messageSlice);
    }
    for (size_t i = 0; i < message.size(); i++) {
        message[i] = abs(message[i]);
    } 

    invertBits(message);
    binaryToString(message);
}

/**
 * Function to decode given input using given parity matrix
*/
void decode(vector<vector<int>> H, vector<int> y) {

    auto m = H.size();
    auto n = H[0].size();

    // Get bits and nodes of the parity-check H matrix
    auto statsVector = bitsNodes(H);
    auto bitsHist = statsVector[0]; auto bits = statsVector[1]; auto nodesHist = statsVector[2]; auto nodes = statsVector[3];

    auto solver = log_belief_propagation;

    auto Lq = zeroMatrix(m, n);
    auto Lr = zeroMatrix(m, n);
    auto Lc = vectorMultiply(y, 2);

    vector<int> x;

    for (int i = 0; i < ITERATION_LIMIT; i++) {
        auto probs = solver(bitsHist, bits, nodesHist, nodes, Lq, Lr, Lc, i);
        Lq = probs.first.first;
        Lr = probs.first.second;
        auto L_post = probs.second;   

        x = {};
        for (auto value : L_post) {
            x.push_back(int(value <= 0));
        }

        auto product = checkMatrixZero(H, x);
        if (product) {
            break;
        }
    }
    
    originalMessage(H, x);
}

void writeCsv(vector<vector<int>> matrix) {
    ofstream f(MATRIX_FILE);
    if (f.is_open()) {
        auto line_size = matrix[0].size();
        for (auto line : matrix) {
            for (size_t i = 0; i < line_size; i++) {
                // If it's the last number, don't put another comma behind it
                if (i == line_size-1) {
                    f << to_string(line[i]);
                }
                else {
                    f << to_string(line[i]) << ",";
                }
            }
            f << '\n';
        }
        f.close();
    }
    else {
        printError("Unable to write the parity check matrix!");
    }
}

vector<vector<int>> readCsv(string filename) {
    vector<vector<int>> matrix;
    string line;
    ifstream f(filename);
    if (f.is_open()) {
        while (getline(f, line)) {
            // Go through characters in the line and numbers between commas wil lbe pushed into the matrix
            string buffer = "";
            vector<int> matrix_line;
            for (auto a : line) {
                if (a == ',') {
                    matrix_line.push_back(stoi(buffer));
                    buffer = "";
                }
                else if (isOneOrZero(a)) {
                    buffer += a;
                }
                // The matrix can only contain 0 or 1
                else {
                    printError("Invalid character in the parity matrix!");
                }
            }
            matrix_line.push_back(stoi(buffer)); // Don't forget the final number which has no comma after
            matrix.push_back(matrix_line);
        }
        f.close();
    }
    else {
        printError("Unable to open the matrix file!");
    }
    return matrix;
}

int main(int argc, char **argv) {
    pair<string, string> args = parseArgs(argc, argv);
    string encode_decode = args.first;
    string matrix_file = args.second;

    vector<vector<int>> user_matrix;

    if (matrix_file != "") {
        user_matrix = readCsv(matrix_file);
    }
    string input = readInput(encode_decode);
    vector<int> binary_input;

    if (encode_decode == "-e") {
        binary_input = encodeBinary(input);
    }
    if (encode_decode == "-d") {
        binary_input = decodeBinary(input);
    }
    if (binary_input.size() == 0) {
        printError("Invalid input! Please enter allowed character only!");
    }

    if (encode_decode == "-e") {
        // If i wasnt given a coding matrix, use my own and write it into the matica.csv
        if (matrix_file == "") {
            pair<vector<vector<int>>, vector<vector<int>>> H_G = make_ldpc(binary_input);
            writeCsv(H_G.first); // Write out the parity matrix (H)
            // Encode using the Generator matrix (G)
            auto encoded_text = encode(H_G.second, binary_input);
            printVector(encoded_text, false);
        }
        // Else use the given one
        else {
            vector<vector<int>> G = construct_G(user_matrix);
            auto encoded_text = encode(G, binary_input);
            printVector(encoded_text, false);
        }
    }
    else if (encode_decode == "-d") {
        // Decode using the given Parity matrix
        decode(user_matrix, binary_input);
    }

    return END_SUCCESS;
}
