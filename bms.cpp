#include <iostream>
#include <string>
#include <bitset> // bitset to convert to binary
#include <vector>
#include <utility> // pair
#include <cmath> // floor
#include <algorithm> // random shuffle

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
        for (auto number : bitset<8>(ascii[i]).to_string()) {
            numbers.push_back(number-48); // ascii 0 is 48, so 48-48 is 0 a 49-48 is 1, just as I want
        }
    }
    return numbers;
}

void printVector(vector<int> vec) {
    for (auto number : vec) {
        cout << number;
    }
    cout << endl;
}

void printVector(vector<vector<int>> vec) {
    for(size_t i = 0; i < vec.size(); i++) {
        for (size_t j = 0; j < vec[0].size(); j++) {
            cout << vec[i][j];
        }   
        cout << endl;
    }
}

vector<vector<int>> zero_matrix(size_t rows, size_t cols) {
    vector<vector<int>> block;
    for (size_t i = 0; i < rows; i++) {
        vector<int> row;
        for (size_t j = 0; j < cols; j++) {
            row.push_back(0);
        }
        block.push_back(row);
    }
    return block;
}

/* Function to simulate python-esque arr[X:Y] */
void setFirstPartOfVector(vector<vector<int>>& arr, int X, int Y, vector<vector<int>> block) {
    auto ctr = 0;
    for (auto i = X; i < Y; i++) {
        arr[i] = block[ctr++];
    }
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

    vector<vector<int>> block = zero_matrix(rows, cols);
    vector<vector<int>> H = zero_matrix(n_equations, codeword);

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

vector<vector<int>> construct_G(vector<vector<int>> H) {
    vector<vector<int>> G;
    return G;
}

/* Function to create LDPC encoding and decoding matrices */
pair<vector<vector<int>>, vector<vector<int>>> make_ldpc(size_t codeword, size_t d_c, size_t d_v) {

    vector<vector<int>> H = construct_H(codeword, d_c, d_v);
    vector<vector<int>> G = construct_G(H);
    return make_pair(H, G);
}

/* Function to encode given string into LDPC code */
void encode(string input) {
    vector<int> encoded = encodeBinary(input);

    size_t codeword = encoded.size() * 2;
    size_t d_c = encoded.size();
    size_t d_v = encoded.size() - 1;

    if (d_v <= 1) {
        printError("Input must be at least 1 character long!");
    }

    pair<vector<vector<int>>, vector<vector<int>>> h_g_pair = make_ldpc(codeword, d_c, d_v);
    
}

void decode(string input) {};

int main(int argc, char **argv) {
    string arg = parseArgs(argc, argv);
    string input = readInput(arg);
    if (arg == "-e") {
        encode(input);
    }
    else if (arg == "-d") {
        decode(input);
    }
    return END_SUCCESS;
}
