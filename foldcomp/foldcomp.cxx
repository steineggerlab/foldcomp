#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <iosfwd>

#include "foldcomp.h"

static PyObject *FoldcompError;

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

int decompress(const char* input, size_t input_size, bool use_alt_order, std::ostream& oss) {
    OneShotReadBuf buf((char*)input, input_size);
    std::istream istr(&buf);

    std::cout.setstate(std::ios_base::failbit);
    Foldcomp compRes;
    int flag = compRes.read(istr);
    std::cout.clear();
    if (flag != 0) {
        return 1;
    }
    std::vector<AtomCoordinate> atomCoordinates;
    compRes.useAltAtomOrder = use_alt_order;
    flag = compRes.decompress(atomCoordinates);
    if (flag != 0) {
        return 1;
    }
    // Write decompressed data to file
    writeAtomCoordinatesToPDB(atomCoordinates, compRes.strTitle, oss);

    return 0;
}


static PyObject *foldcomp_decompress(PyObject *self, PyObject *args) {
    // Unpack a string from the arguments
    const char *strArg;
    Py_ssize_t strSize;
    if (!PyArg_ParseTuple(args, "y#", &strArg, &strSize)) {
        return NULL;
    }

    std::ostringstream oss;
    int err = decompress(strArg, strSize, true, oss);
    if (err != 0) {
        PyErr_SetString(FoldcompError, "Error decompressing");
        return NULL;
    }

    PySys_WriteStdout("%s %s\n", strArg, oss.str().c_str());

    return PyBytes_FromStringAndSize(oss.str().c_str(), oss.str().size());
}

static PyMethodDef foldcomp_methods[] = {
    // {"compress", foldcomp_compress, METH_VARARGS, "Compress a PDB file."},
    {"decompress", foldcomp_decompress, METH_VARARGS, "Decompress a PDB file."},

    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef foldcomp_module_def = {
    PyModuleDef_HEAD_INIT,
    "foldcomp",
    NULL,
    -1,
    foldcomp_methods
};

PyMODINIT_FUNC PyInit_foldcomp(void) {
  PyObject *m = PyModule_Create(&foldcomp_module_def);
    if (m == NULL) {
        return NULL;
    }

    FoldcompError = PyErr_NewException("foldcomp.error", NULL, NULL);
    Py_XINCREF(FoldcompError);
    if (PyModule_AddObject(m, "error", FoldcompError) < 0) {
        Py_XDECREF(FoldcompError);
        Py_CLEAR(FoldcompError);
        Py_DECREF(m);
        return NULL;
    }

  return m;
}
