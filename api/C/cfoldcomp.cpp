#include "foldcomp.h"
#include <iosfwd>

int compress(const char* buffer, char* output, int* output_size) {
    return 0;
}

// https://stackoverflow.com/questions/1448467/initializing-a-c-stdistringstream-from-an-in-memory-buffer/1449527
struct OneShotReadBuf : public std::streambuf
{
    OneShotReadBuf(char* s, std::size_t n)
    {
        setg(s, s, s + n);
    }
};

int decompress(const char* input, size_t input_size, bool use_alt_order, char* output, int* output_size) {
    OneShotReadBuf buf((char*)input, input_size);
    std::istream istr(&buf);

    int flag = 0;
    Foldcomp compRes = Foldcomp();
    flag = compRes.read(istr);
    if (flag != 0) {
        std::cerr << "Error reading" << std::endl;
        return 1;
    }
    std::vector<AtomCoordinate> atomCoordinates;
    compRes.useAltAtomOrder = use_alt_order;
    flag = compRes.decompress(atomCoordinates);
    if (flag != 0) {
        std::cerr << "Error decompressing compressed data." << std::endl;
        return 1;
    }
    // Write decompressed data to file
    std::ostringstream oss;
    writeAtomCoordinatesToPDB(atomCoordinates, compRes.strTitle, oss);

    std::string pdb = oss.str();
    output = (char*) malloc(pdb.size() * sizeof(char));
    if (output == NULL) {
        return 0;
    }
    std::copy(pdb.begin(), pdb.end(), output);

    return 0;
}