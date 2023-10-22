#include <iostream>
#include <string>

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

void encode(string input) {
    cout << input << endl;
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
