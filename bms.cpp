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

/* Function to modulo each value in the matrix by a given number */
void moduloMatrix(vector<vector<int>> &arr, int n) {
    for(size_t i = 0; i < arr.size(); i++) {
        for(size_t j = 0; j < arr[i].size(); j++) {
            arr[i][j] = arr[i][j] % n;
        }
    }
}

/* Function to modulo each value in the vector by a given number */
void moduloMatrix(vector<int> &arr, int n) {
    for(size_t i = 0; i < arr.size(); i++) {
            arr[i] = (n + (arr[i] % n)) % n;
    }
}

/* Function to calculate binary product of 2 matrices */
vector<int> binaryProduct(vector<vector<int>> H, vector<int> x) {
    auto prod = matrixProduct(H, x);
    moduloMatrix(prod, 2);
    return prod;
}

/** Function to calculate binary product of 2 matrices */
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

/** Function to print the vector to stdout */
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

/** Function to print the vector to stdout */
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
        // zero padding makes it unable to easily get ascii values back - dont use it
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


/* Function to print information about the vector */
void vectorInfo(vector<int> vec) {
    cout << vec.size() << endl;
}

/* Function to print information about the vector */
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

/* Function to get python-like slice of given matrix */
vector<int> getSlice(vector<vector<int>> arr, int X, int Y, int Z) {
    vector<int> slice;
    for (size_t i = X; (int) i < Y; i++) {
        slice.push_back(arr[i][Z]);
    }
    return slice;
}

/* Function to get python-like slice of given matrix */
vector<int> getSlice(vector<int> arr, int X, int Y) {
    vector<int> slice;
    for (size_t i = X; (int) i < Y; i++) {
        slice.push_back(arr[i]);
    }
    return slice;
}

/* Fuction to transpose a matrix in vector form */
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
vector<vector<int>> createParityMatrix(size_t codeword, size_t d_c, size_t d_v) {

    auto eq_count = floor((codeword * d_v) / d_c);
    auto bs = floor(eq_count / d_v);

    size_t rows = floor(eq_count / d_v);
    size_t cols = codeword;

    vector<vector<int>> block = zeroMatrix(rows, cols);
    vector<vector<int>> parityMatrix = zeroMatrix(eq_count, codeword);

    // Fill the first block with consecutive ones in each block row
    for (auto i = 0; i < bs; i++) {
        for (auto j = i * d_c; j < (i+1) * d_c; j++) {
            block[i][j] = 1;
        }
    }
    // set the H vector up to the block size to correspond with the set block
    setFirstPartOfVector(parityMatrix, 0, bs, block);

    // remaining blocks are permutations of the first block's column
    for (size_t i = 1; i < d_v; i++) {
        transpose(block);
        random_shuffle(block.begin(), block.end());
        transpose(block);
        setFirstPartOfVector(parityMatrix, i * bs, (i+1) * bs, block);
    }
    return parityMatrix;   
}

/* Function to simulate argMax() in numpy */
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

/* Function to return identity matrix of given size */
vector<vector<int>> identityVector(size_t n) {
    vector<vector<int>> identity = zeroMatrix(n, n);
    for (size_t i = 0; i < n; i++) {
        identity[i][i] = 1;
    }
    return identity;
}

/* Function to calculate gauss elimination */
pair<vector<vector<int>>, vector<vector<int>>> gaussjordan(vector<vector<int>> matrix, bool toggle) {
    
    // initial matrix size
    auto rows = matrix.size(); // rows
    auto cols = matrix[0].size(); // cols

    size_t prev_piv = -1;

    vector<vector<int>> p_matrix;
    
    if (toggle) {
        p_matrix = identityVector(rows);
    }

    // go through each column
    for (size_t j = 0; j < cols; j++) {
        vector<int> filtre_down = getSlice(matrix, prev_piv+1, (int) rows, (int) j); // part of current column below the old pivot
        size_t pivot = largestvalueIndex(filtre_down) + prev_piv + 1; // new pivot becomes the largest value index in the selected part of the column

        // if pivot is 1
        if (matrix[pivot][j]) {
            prev_piv++;

            // swap rows if pivot index has moved by more than one
            if (prev_piv != pivot) {
                auto aux = matrix[pivot];
                matrix[pivot] = matrix[prev_piv];
                matrix[prev_piv] = aux;

                if (toggle) {
                    aux = p_matrix[pivot];
                    p_matrix[pivot] = p_matrix[prev_piv];
                    p_matrix[prev_piv] = aux;
                }
            }
            
            // remove other values in current column
            for (size_t i = 0; i < rows; i++) {
                if (i != prev_piv && matrix[i][j]) {
                    if (toggle) {
                        vector<int> tmp = {};
                        for (size_t q = 0; q < p_matrix[0].size(); q++) {
                            auto abs_val = abs(p_matrix[i][q] - p_matrix[prev_piv][q]);
                            tmp.push_back(abs_val);
                        }
                        p_matrix[i] = tmp;
                    }
                    vector<int> tmp = {};
                        for (size_t q = 0; q < matrix[0].size(); q++) {
                            auto abs_val = abs(matrix[i][q] - matrix[prev_piv][q]);
                            tmp.push_back(abs_val);
                        }
                    matrix[i] = tmp;
                }
            }
        }
        // break after going through each row
        if (prev_piv == rows-1) {
            break;
        }
    }
    return make_pair(matrix, p_matrix);
}

/* Function to calculate the sum of a matrix, numpy-like sum() */
size_t vecSum(vector<int> arr) {
    size_t sum = 0;
    for (size_t k = 0; k < arr.size(); k++) {
        sum += arr[k];
    }
    return sum;
}

/* Function to calculate the sum of a matrix, numpy-like sum() */
size_t vecSum(vector<vector<int>> arr) {
    size_t sum = 0;
    for (size_t i = 0; i < arr.size(); i++) {
        for (size_t k = 0; k < arr[i].size(); k++) {
            sum += arr[i][k];
        }
    }
    return sum;
}

/* Function to multiply all values in a vector by a constant */
vector<double> vectorMultiply(vector<double> &A, double val) {
    vector<double> multiplied;
    for (size_t i = 0; i < A.size(); i++) {
        multiplied.push_back(A[i] * val);
    } 
    return multiplied;
}

/* Function to multiply all values in a vector by a constant */
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

/* Function to construct the generating matrix from the parity matrix */
vector<vector<int>> createGeneratingMatrix(vector<vector<int>> parityMatrix) {

    auto codewordSize = parityMatrix[0].size(); // codeword size (n of columns in the parity matrix)

    transpose(parityMatrix);
    // do gauss-jordan elimination on transposed H
    auto H_tQ = gaussjordan(parityMatrix, true);
    transpose(parityMatrix);

    auto parityReduced = H_tQ.first; // parity matrix in reduced row echelon form

    transpose(parityReduced);

    // do gaussjordan on the transposed reduced row echelon form of parity matrix
    auto diag = gaussjordan(parityReduced, false).first; // diagonalized form of the reduced row echelon form
    transpose(parityReduced);

    auto transformationMatrix = H_tQ.second;
    transpose(transformationMatrix);

    auto bitCount = codewordSize - vecSum(diag); // number of information bits in the codeword

    auto tmp = zeroMatrix(codewordSize, bitCount);

    setFirstPartOfVector(tmp, codewordSize - bitCount, (int) tmp.size(), identityVector(bitCount));

    auto generatingMatrix = matrixProduct(transformationMatrix, tmp);
    moduloMatrix(generatingMatrix, 2);

    return generatingMatrix;
}

/* Function to create LDPC Parity matrix and generating matrix */
pair<vector<vector<int>>, vector<vector<int>>> createMatrixes(vector<int> binary_input) {

    size_t codeword = binary_input.size() * 2;
    size_t d_c = binary_input.size();
    size_t d_v = binary_input.size() - 1;

    if (d_v <= 1) {
        printError("Input must be at least 1 character long!");
    }

    vector<vector<int>> parity = createParityMatrix(codeword, d_c, d_v);
    vector<vector<int>> generating = createGeneratingMatrix(parity);

    return make_pair(parity, generating);
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


/* Function to encode the binary string by multiplying it with the generating matrix */
vector<int> encode(vector<vector<int>> generatingMatrix, vector<int> message) {
    auto encoded = binaryProduct(generatingMatrix, message);
    return encoded;
}

/* Function to simulate numpy where */
vector<vector<int>> where(vector<vector<int>> matrix) {
    vector<int> row;
    vector<int> col;
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t k = 0; k < matrix[i].size(); k++) {
            if(matrix[i][k]) {
                col.push_back(k);
                row.push_back(i);
            }
        }
    }
    vector<vector<int>> rowCols = {row, col};
    return rowCols;
}

/* Function to count number of occurences of a given number in a given array */
int count(vector<int> arr, int x) {
    int ctr = 0;
    for (auto number : arr) {
        if (number == x) {
            ctr++;
        }
    }
    return ctr;
}

/* Function to simulate the binCount numpy function */
vector<int> binCount(vector<int> arr) {
    vector<int> bins;

    int max = *max_element(arr.begin(), arr.end());

    for (int i = 0; i < max+1; i++) {
        auto n_of_occurences = count(arr, i);
        bins.push_back(n_of_occurences);
    }
    return bins;
}

/* Function to find the bits and nodes of the parity matrix */
vector<vector<int>> bitsNodes(vector<vector<int>> parityMatrix) {

    auto bitIndexes = where(parityMatrix); // indexes of 1s
    auto bitRows = bitIndexes[0];
    auto bitCols = bitIndexes[1];

    transpose(parityMatrix);

    auto nodeIndexes = where(parityMatrix); // indexes of 1s of transposed parity (nodes)
    auto nodeRows = nodeIndexes[0];
    auto nodeCols = nodeIndexes[1];

    auto bitsBins = binCount(bitRows); // occurences of bits
    auto nodeBins = binCount(nodeRows); // occurences of nodes

    vector<vector<int>> bitsNodes = {bitsBins, bitCols, nodeBins, nodeCols};
    return bitsNodes;
}

/* Function to calculate the probabilities using belief propagation 
 * aprioriMatrix is the apriori information
 * CVM is check:variable matrix
 * VCM is variable:check matrix
*/
pair<pair<vector<vector<int>>, vector<vector<int>>>, vector<double>> log_belief_propagation(vector<int> bitsBins, vector<int> bitCols, 
                                                                          vector<int> nodesBins, vector<int> nodeCols, 
                                                                          vector<vector<int>> VCM, vector<vector<int>> CVM, vector<int> aprioriMatrix, int toggle) {

    auto rows = CVM.size(); // rows
    auto cols = CVM[0].size(); // cols
    auto bits_counter = 0;
    auto nodes_counter = 0;

    // Horizontal
    for (size_t row = 0; row < rows; row++) {
        auto currBin = bitsBins[row];
        auto currBitColSlice = getSlice(bitCols, bits_counter, bits_counter + currBin);
        bits_counter += currBin;
        for (auto bitCol : currBitColSlice) {
            auto tmpCurrBitColSlice = currBitColSlice;
            int val = 1;
            // if first iterationm multiply by apriori information
            if (toggle == 0) {
                for (size_t k = 0; k < tmpCurrBitColSlice.size(); k++) {
                    if (tmpCurrBitColSlice[k] != bitCol) {
                        val *= tanh(0.5 * aprioriMatrix[tmpCurrBitColSlice[0]]);
                    }
                }
            }
            // else multiply by VCM
            else {
                for (size_t k = 0; k < tmpCurrBitColSlice.size(); k++) {
                    if (tmpCurrBitColSlice[k] != bitCol) {
                        val *= tanh(0.5 * VCM[row][tmpCurrBitColSlice[0]]);
                    }
                }
            }
            auto numerator = val + 1;
            auto denominator = 1 - val;
            if (numerator == 0) {
                CVM[row][bitCol] = -1;
            }
            else if (denominator == 0) {
                CVM[row][bitCol] = 1;
            }
            else {
                CVM[row][bitCol] = log(numerator / denominator);
            }
        }
    }

    // Vertical
    for (size_t col = 0; col < cols; col++) {
        auto currBin = nodesBins[col];
        auto currNodeColSlice = getSlice(nodeCols, nodes_counter, nodes_counter + currBin);
        nodes_counter += currBin;
        for (auto currCol : currNodeColSlice) {
            auto tmpCurrNodeColSlice = currNodeColSlice;
            VCM[currCol][col] = aprioriMatrix[col];
            for (size_t k = 0; k < tmpCurrNodeColSlice.size(); k++) {
                if (tmpCurrNodeColSlice[k] != currCol) {
                    VCM[currCol][col] += CVM[tmpCurrNodeColSlice[k]][col];
                }
            }
        }
    }

    // Posterior
    vector<double> posteriorMatrix;
    for (size_t col = 0; col < cols; col++) {
        posteriorMatrix.push_back(0);
    }
    nodes_counter = 0;

    for (size_t col = 0; col < cols; col++) {
        auto nodeBin = nodesBins[col];
        auto currNodeColSlice = getSlice(nodeCols, nodes_counter, nodes_counter + nodeBin);
        nodes_counter += nodeBin;
        int value = 0;
        for (auto number : currNodeColSlice) {
            value += CVM[number][col];
        }
        posteriorMatrix[col] = aprioriMatrix[col] + value;
    }

    return make_pair(make_pair(VCM, CVM), posteriorMatrix);
}

/* Function to check if binary product of given matrix contains only zeros */
int checkMatrixZero(vector<vector<int>> matrix, vector<int> vec) {
    auto prod = binaryProduct(matrix, vec);
    for (auto a : prod) {
        if (a != 0) {
            return 0;
        }
    }
    return 1;
}

/* Function to perform gauss elimination */
pair<vector<vector<int>>, vector<int>> gaussElimination(vector<vector<int>> matrix, vector<int> vec) {
    auto rows = matrix.size(); // rows
    auto cols = matrix[0].size(); // cols

    auto limit = rows < cols ? rows : cols; // lower number is the limit

    // go through columns
    for (size_t col = 0; col < limit; col++) {
        vector<int> pivots;

        // find row indexes with 1s in current column
        for (size_t i = col; i < rows; i++) {
            if (matrix[i][col]) {
                pivots.push_back(i);
            }
        }

        size_t pivot = 99999999;

        // If there was at least one 1 found in the column, base the pivot on the lowest value
        if (pivots.size()) {
            pivot = *min_element(pivots.begin(), pivots.end());
        }
        else {
            continue;
        }

        // if pivot doesnt correspon to current column, swap them
        if (pivot != col) {
            auto tmp = matrix[col];
            matrix[col] = matrix[pivot];
            matrix[pivot] = tmp;

            auto tmp2 = vec[col];
            vec[col] = vec[pivot];
            vec[pivot] = tmp2;
        }

        // do the elimination
        for (auto i = col+1; i < rows; i++) {
            if (matrix[i][col]) {
                // substract rows from following rows to remove 1s
                for (size_t k = 0; k < matrix[0].size(); k++) {
                    matrix[i][k] = abs(matrix[i][k] - matrix[col][k]);
                }
                vec[i] = abs(vec[i] - vec[col]);
            }
        }
    }
    return make_pair(matrix, vec);
}

/* Function to invert bits in given vector */
void invertBits(vector<int>&vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        vec[i] = vec[i] == 0 ? 1 : 0;
    }
}

/* Function to print binary string representation as ascii string*/
void binaryToString(vector<int> &vec) {

    string tmp = "";

    for (size_t i = 0; i < vec.size(); i++) {
        if (i % 8 == 0 && i > 0) {
            auto ascii_val = stoi(tmp, nullptr, 2);
            cout << (char) ascii_val;
            tmp = to_string(vec[i]);
        }
        else {
            tmp += to_string(vec[i]);
        }
    }
    auto ascii_val = stoi(tmp, nullptr, 2);
    cout << (char) ascii_val; // Finish the last one
    cout << endl;
}

/* Function to get the original message from the decoded sequence and print it */
void originalMessage(vector<vector<int>> parityMatrix, vector<int> msg) {
    vector<vector<int>> generatingMatrix = createGeneratingMatrix(parityMatrix); // find generating matrix from parity

    auto generatingCols = generatingMatrix[0].size(); // information bits in original msg
    auto gaussResult = gaussElimination(generatingMatrix, msg);
    auto eliminatedMatrix = gaussResult.first;
    auto eliminatedMsg = gaussResult.second;

    vector<int> message;
    for (size_t i = 0; i < generatingCols; i++) {
        message.push_back(0);
    }

    message[generatingCols-1] = eliminatedMsg[generatingCols-1];
    for (int i = generatingCols - 2; i >= 0; i--) {
        message[i] = eliminatedMsg[i];

        vector<int> eliminatedSlice;
        vector<int> messageSlice;
        for (size_t z = i+1; z < generatingCols; z++) {
            eliminatedSlice.push_back(eliminatedMatrix[i][z]);
            messageSlice.push_back(message[z]);
        }
        message[i] -= binaryProduct(eliminatedSlice, messageSlice);
    }
    // all elements must be positive
    for (size_t i = 0; i < message.size(); i++) {
        message[i] = abs(message[i]);
    } 

    invertBits(message);

    // convert bits to string and print them
    binaryToString(message);
}

/**
 * Function to decode given input using given parity matrix
*/
void decode(vector<vector<int>> parityMatrix, vector<int> received_msg) {

    auto rows = parityMatrix.size();
    auto cols = parityMatrix[0].size();

    // Get bits and nodes of the parity-check H matrix
    auto statsVector = bitsNodes(parityMatrix);
    auto bitsBins = statsVector[0]; auto bitCols = statsVector[1]; auto nodesBins = statsVector[2]; auto nodeCols = statsVector[3];
    auto VCM = zeroMatrix(rows, cols);
    auto CVM = zeroMatrix(rows, cols);
    auto apriori = received_msg; // given vector is the apriori information

    vector<int> orig_msg;

    for (int i = 0; i < ITERATION_LIMIT; i++) {
        auto probs = log_belief_propagation(bitsBins, bitCols, nodesBins, nodeCols, VCM, CVM, apriori, i);
        VCM = probs.first.first;
        CVM = probs.first.second;
        auto posterior = probs.second;   

        orig_msg = {};
        for (auto value : posterior) {
            orig_msg.push_back(int(value <= 0));
        }

        auto check = checkMatrixZero(parityMatrix, orig_msg);
        if (check) { // early convergence
            break;
        }
    }
    
    originalMessage(parityMatrix, orig_msg); // get the original message from decoded codeword
}

/* Function to write a given matrix as a csv file */
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

/* Function to load parity matrix from given csv file */
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

/* Main */
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
            pair<vector<vector<int>>, vector<vector<int>>> H_G = createMatrixes(binary_input);
            writeCsv(H_G.first); // Write out the parity matrix (H)
            // Encode using the Generator matrix (G)
            auto encoded_text = encode(H_G.second, binary_input);
            printVector(encoded_text, false);
        }
        // Else use the given one
        else {
            vector<vector<int>> G = createGeneratingMatrix(user_matrix);
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
