#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <sstream> // IWYU pragma: keep

#include "atom_coordinate.h"
#include "foldcomp.h"
#include "database_reader.h"

static PyObject *FoldcompError;

typedef struct {
    PyObject_HEAD
    PyObject* user_ids;
    bool decompress;
    void* memory_handle;
} FoldcompDatabaseObject;

int decompress(const char* input, size_t input_size, bool use_alt_order, std::ostream& oss, std::string& name);
static PyObject* FoldcompDatabase_close(PyObject* self);
static PyObject* FoldcompDatabase_enter(PyObject* self);
static PyObject* FoldcompDatabase_exit(PyObject* self, PyObject* args);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wcast-function-type"
static PyMethodDef FoldcompDatabase_methods[] = {
    {"close", (PyCFunction)FoldcompDatabase_close, METH_NOARGS, "Close the database."},
    {"__enter__", (PyCFunction)FoldcompDatabase_enter, METH_NOARGS, "Enter the runtime context related to this object."},
    {"__exit__", (PyCFunction)FoldcompDatabase_exit, METH_VARARGS, "Exit the runtime context related to this object."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};
#pragma GCC diagnostic pop

// FoldcompDatabase_sq_length
static Py_ssize_t FoldcompDatabase_sq_length(PyObject* self) {
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;
    if (db->user_ids != NULL) {
        return PySequence_Length(db->user_ids);
    }
    return (Py_ssize_t)reader_get_size(db->memory_handle);
}

// FoldcompDatabase_sq_item
static PyObject* FoldcompDatabase_sq_item(PyObject* self, Py_ssize_t index) {
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;

    const char* data;
    size_t length;
    if (db->user_ids != NULL) {
        if (index >= PySequence_Length(db->user_ids)) {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }
        PyObject* item = PySequence_GetItem(db->user_ids, index);
        // get string representation of id as c string
        uint32_t key = reader_lookup_entry(db->memory_handle, PyUnicode_AsUTF8(item));
        if (key == UINT32_MAX) {
            PyErr_SetString(PyExc_KeyError, "Could not find key in database.");
            return NULL;
        }
        int64_t id = reader_get_id(db->memory_handle, key);
        if (id == -1) {
            PyErr_SetString(PyExc_KeyError, "Could not find key in database.");
            return NULL;
        }
        data = reader_get_data(db->memory_handle, id);
        length = std::max(reader_get_length(db->memory_handle, id), (int64_t)1) - (int64_t)1;
    } else {
        if (index >= reader_get_size(db->memory_handle)) {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }
        data = reader_get_data(db->memory_handle, index);
        length = std::max(reader_get_length(db->memory_handle, index), (int64_t)1) - (int64_t)1;
    }
    if (db->decompress) {
        std::ostringstream oss;
        std::string name;
        int err = decompress(data, length, false, oss, name);
        if (err != 0) {
            PyErr_SetString(FoldcompError, "Error decompressing.");
            return NULL;
        }
        return Py_BuildValue("(s,O)", name.c_str(), PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, oss.str().c_str(), oss.str().size()));
    }
    return PyBytes_FromStringAndSize(data, length);
}

// PySequenceMethods
static PySequenceMethods FoldcompDatabase_as_sequence = {
    &FoldcompDatabase_sq_length, // sq_length
    0, // sq_concat
    0, // sq_repeat
    &FoldcompDatabase_sq_item, // sq_item
    0, // sq_slice
    0, // sq_ass_item
    0, // sq_ass_slice
    0, // sq_contains
    0, // sq_inplace_concat
    0, // sq_inplace_repeat
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
static PyTypeObject FoldcompDatabaseType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "foldcomp.FoldcompDatabase",    /* tp_name */
    sizeof(FoldcompDatabaseObject), /* tp_basicsize */
    0,                         /* tp_itemsize */
    0,                         /* tp_dealloc */
    0,                         /* tp_vectorcall_offset */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_as_async */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    &FoldcompDatabase_as_sequence, /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "FoldcompDatabase objects", /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    FoldcompDatabase_methods,  /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    0,                         /* tp_new */
    0,                         /* tp_free */
    0,                         /* tp_is_gc */
    0,                         /* tp_bases */
    0,                         /* tp_mro */
    0,                         /* tp_cache */
    0,                         /* tp_subclasses */
    0,                         /* tp_weaklist */
    0,                         /* tp_del */
    0,                         /* tp_version_tag */
    0,                         /* tp_finalize */
    //0,                         /* tp_vectorcall */
};
#pragma GCC diagnostic pop

// FoldcompDatabase_close
static PyObject* FoldcompDatabase_close(PyObject* self) {
    if (!PyObject_TypeCheck(self, &FoldcompDatabaseType)) {
        PyErr_SetString(PyExc_TypeError, "Expected FoldcompDatabase object.");
        return NULL;
    }
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;
    Py_XDECREF(db->user_ids);
    if (db->memory_handle != NULL) {
        free_reader(db->memory_handle);
        db->memory_handle = NULL;
    }
    Py_RETURN_NONE;
}

// FoldcompDatabase_enter
static PyObject* FoldcompDatabase_enter(PyObject* self) {
    Py_INCREF(self);
    return (PyObject*)self;
}

// FoldcompDatabase_exit
static PyObject *FoldcompDatabase_exit(PyObject *self, PyObject* /* args */) {
    return FoldcompDatabase_close(self);
}

// https://stackoverflow.com/questions/1448467/initializing-a-c-stdistringstream-from-an-in-memory-buffer/1449527
struct OneShotReadBuf : public std::streambuf
{
    OneShotReadBuf(char* s, std::size_t n)
    {
        setg(s, s, s + n);
    }
};

// Decompress
int decompress(const char* input, size_t input_size, bool use_alt_order, std::ostream& oss, std::string& name) {
    OneShotReadBuf buf((char*)input, input_size);
    std::istream istr(&buf);

    std::cout.setstate(std::ios_base::failbit);
    Foldcomp compRes;
    int flag = compRes.read(istr);
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
    std::cout.clear();

    name = compRes.strTitle;

    return 0;
}
// Python binding for decompress
static PyObject *foldcomp_decompress(PyObject* /* self */, PyObject *args) {
    // Unpack a string from the arguments
    const char *strArg;
    Py_ssize_t strSize;
    if (!PyArg_ParseTuple(args, "y#", &strArg, &strSize)) {
        return NULL;
    }

    std::ostringstream oss;
    std::string name;
    int err = decompress(strArg, strSize, false, oss, name);
    if (err != 0) {
        PyErr_SetString(FoldcompError, "Error decompressing.");
        return NULL;
    }

    return Py_BuildValue("(s,O)", name.c_str(), PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, oss.str().c_str(), oss.str().size()));
}

std::string trim(const std::string& str, const std::string& whitespace = " \t") {
    const std::string::size_type strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const std::string::size_type strEnd = str.find_last_not_of(whitespace);
    const std::string::size_type strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

// Compress
int compress(const std::string& name, const std::string& pdb_input, std::ostream& oss, int anchor_residue_threshold) {
    std::vector<AtomCoordinate> atomCoordinates;
    // parse ATOM lines from PDB file into atomCoordinates
    std::istringstream iss(pdb_input);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.substr(0, 4) == "ATOM") {
            atomCoordinates.emplace_back(
                trim(line.substr(12, 4)), // atom
                trim(line.substr(17, 3)), // residue
                line.substr(21, 1), // chain
                std::stoi(line.substr(6,  5)), // atom_index
                std::stoi(line.substr(22, 4)), // residue_index
                std::stof(line.substr(30, 8)), std::stof(line.substr(38, 8)), std::stof(line.substr(46, 8)), // coordinates
                std::stof(line.substr(54, 6)), // occupancy
                std::stof(line.substr(60, 6)) // tempFactor
            );
        }
    }
    if (atomCoordinates.size() == 0) {
        return 1;
    }

    removeAlternativePosition(atomCoordinates);

    // compress
    Foldcomp compRes;
    compRes.strTitle = name;
    compRes.anchorThreshold = anchor_residue_threshold;
    compRes.compress(atomCoordinates);
    compRes.writeStream(oss);

    return 0;
}
// Python binding for compress
static PyObject *foldcomp_compress(PyObject* /* self */, PyObject *args, PyObject* kwargs) {
    const char* name;
    const char* pdb_input;
    PyObject* anchor_residue_threshold = NULL;
    static const char *kwlist[] = {"name", "pdb_content", "anchor_residue_threshold", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|$O", const_cast<char**>(kwlist), &name, &pdb_input, &anchor_residue_threshold)) {
        return NULL;
    }

    if (anchor_residue_threshold != NULL && !PyLong_Check(anchor_residue_threshold)) {
        PyErr_SetString(PyExc_TypeError, "anchor_residue_threshold must be an integer");
        return NULL;
    }

    int threshold = DEFAULT_ANCHOR_THRESHOLD;
    if (anchor_residue_threshold != NULL) {
        threshold = PyLong_AsLong(anchor_residue_threshold);
    }

    std::ostringstream oss;
    int flag = compress(name, pdb_input, oss, threshold);
    if (flag != 0) {
        PyErr_SetString(FoldcompError, "Error compressing.");
        return NULL;
    }

    return PyBytes_FromStringAndSize(oss.str().c_str(), oss.str().length());
}


PyTypeObject* pathType = NULL;

static PyObject *foldcomp_open(PyObject* /* self */, PyObject* args, PyObject* kwargs) {
    PyObject* path;
    PyObject* user_ids = NULL;
    PyObject* decompress = NULL;
    static const char *kwlist[] = {"path", "ids", "decompress", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&|$OO", const_cast<char**>(kwlist), PyUnicode_FSConverter, &path, &user_ids, &decompress)) {
        return NULL;
    }
    if (path == NULL) {
        PyErr_SetString(PyExc_TypeError, "path must be a path-like object");
        return NULL;
    }
    const char* pathCStr = PyBytes_AS_STRING(path);
    if (pathCStr == NULL) {
        Py_XDECREF(path);
        PyErr_SetString(PyExc_TypeError, "path must be a path-like object");
        return NULL;
    }

    if (user_ids != NULL && !PyList_Check(user_ids)) {
        Py_XDECREF(path);
        PyErr_SetString(PyExc_TypeError, "user_ids must be a list.");
        return NULL;
    }

    if (decompress != NULL && !PyBool_Check(decompress)) {
        Py_XDECREF(path);
        PyErr_SetString(PyExc_TypeError, "decompress must be a boolean");
        return NULL;
    }

    std::string dbname(pathCStr);
    std::string index = dbname + ".index";
    Py_XDECREF(path);

    FoldcompDatabaseObject *obj = PyObject_New(FoldcompDatabaseObject, &FoldcompDatabaseType);
    if (obj == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for FoldcompDatabaseObject");
        return NULL;
    }

    obj->user_ids = NULL;
    int mode = DB_READER_USE_DATA;
    if (user_ids != NULL && PySequence_Length(user_ids) > 0) {
        mode |= DB_READER_USE_LOOKUP;
        obj->user_ids = user_ids;
        Py_INCREF(obj->user_ids);
    }

    if (decompress == NULL) {
        obj->decompress = true;
    } else {
        obj->decompress = PyObject_IsTrue(decompress);
    }

    obj->memory_handle = make_reader(dbname.c_str(), index.c_str(), mode);
    return (PyObject*)obj;
}

// C++ vector to Python list
// Original code from https://gist.github.com/rjzak/5681680

PyObject* vectorToList_Float(const std::vector<float>& data) {
    PyObject* listObj = PyList_New(data.size());
    if (!listObj) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for list");
        return NULL;
    }
    for (size_t i = 0; i < data.size(); i++) {
        PyObject* num = PyFloat_FromDouble((double)data[i]);
        if (!num) {
            Py_DECREF(listObj);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for list");
            return NULL;
        }
        PyList_SET_ITEM(listObj, i, num);
    }
    return listObj;
}

PyObject* vector2DToList_Float(const std::vector<float3d>& data) {
    PyObject* listObj = PyList_New(data.size());
    if (!listObj) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for list");
        return NULL;
    }
    for (size_t i = 0; i < data.size(); i++) {
        PyObject* inner = Py_BuildValue("(f,f,f)", data[i].x, data[i].y, data[i].z);
        if (!inner) {
            Py_DECREF(listObj);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for list");
            return NULL;
        }
        PyList_SET_ITEM(listObj, i, inner);
    }
    return listObj;
}

PyObject* getPyDictFromFoldcomp(Foldcomp* fcmp, const std::vector<float3d>& coords) {
    // Output: Dictionary
    PyObject* dict = PyDict_New();
    if (dict == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Python dictionary");
        return NULL;
    }
    PyObject* result = NULL;

    // Dictionary keys: phi, psi, omega, torsion_angles, residues, bond_angles, coordinates
    // Convert vectors to Python lists
    PyObject* phi = vectorToList_Float(fcmp->phi);
    if (phi == NULL) {
        Py_XDECREF(dict);
        return NULL;
    }
    PyObject* psi = vectorToList_Float(fcmp->psi);
    if (psi == NULL) {
        Py_XDECREF(dict);
        Py_XDECREF(phi);
        return NULL;
    }
    PyObject* omega = vectorToList_Float(fcmp->omega);
    if (omega == NULL) {
        Py_XDECREF(dict);
        Py_XDECREF(phi);
        Py_XDECREF(psi);
        return NULL;
    }
    PyObject* torsion_angles = vectorToList_Float(fcmp->backboneTorsionAngles);
    if (torsion_angles == NULL) {
        Py_XDECREF(dict);
        Py_XDECREF(phi);
        Py_XDECREF(psi);
        Py_XDECREF(omega);
        return NULL;
    }
    PyObject* bond_angles = vectorToList_Float(fcmp->backboneBondAngles);
    if (bond_angles == NULL) {
        Py_XDECREF(dict);
        Py_XDECREF(phi);
        Py_XDECREF(psi);
        Py_XDECREF(omega);
        Py_XDECREF(torsion_angles);
        return NULL;
    }
    PyObject* residues = PyUnicode_FromStringAndSize(fcmp->residues.data(), fcmp->residues.size());
    if (residues == NULL) {
        Py_XDECREF(dict);
        Py_XDECREF(phi);
        Py_XDECREF(psi);
        Py_XDECREF(omega);
        Py_XDECREF(torsion_angles);
        Py_XDECREF(bond_angles);
        return NULL;
    }
    PyObject* b_factors = vectorToList_Float(fcmp->tempFactors);
    if (b_factors == NULL) {
        Py_XDECREF(dict);
        Py_XDECREF(phi);
        Py_XDECREF(psi);
        Py_XDECREF(omega);
        Py_XDECREF(torsion_angles);
        Py_XDECREF(bond_angles);
        Py_XDECREF(residues);
        return NULL;
    }

    PyObject* coordinates = vector2DToList_Float(coords);
    if (coordinates == NULL) {
        Py_XDECREF(dict);
        Py_XDECREF(phi);
        Py_XDECREF(psi);
        Py_XDECREF(omega);
        Py_XDECREF(torsion_angles);
        Py_XDECREF(bond_angles);
        Py_XDECREF(residues);
        Py_XDECREF(b_factors);
        return NULL;
    }

    // Set dictionary keys and values
    PyDict_SetItemString(dict, "phi", phi);
    PyDict_SetItemString(dict, "psi", psi);
    PyDict_SetItemString(dict, "omega", omega);
    PyDict_SetItemString(dict, "torsion_angles", torsion_angles);
    PyDict_SetItemString(dict, "bond_angles", bond_angles);
    PyDict_SetItemString(dict, "residues", residues);
    PyDict_SetItemString(dict, "b_factors", b_factors);
    PyDict_SetItemString(dict, "coordinates", coordinates);

    result = dict;

    if (result == NULL) {
        Py_XDECREF(dict);
        Py_XDECREF(phi);
        Py_XDECREF(psi);
        Py_XDECREF(omega);
        Py_XDECREF(torsion_angles);
        Py_XDECREF(bond_angles);
        Py_XDECREF(residues);
        Py_XDECREF(b_factors);
        Py_XDECREF(coordinates);
        return NULL;
    }

    return result;
}

// Extract
// Return a dictionary with the following keys:
// phi, psi, omega, torsion_angles, residues, bond_angles, coordinates, b_factors
// 01. Extract information starting from FCZ file
PyObject* getDataFromFCZ(const char* input, size_t input_size) {
    // Input
    OneShotReadBuf buf((char*)input, input_size);
    std::istream istr(&buf);

    Foldcomp compRes;
    int flag = compRes.read(istr);
    if (flag != 0) {
        PyErr_SetString(PyExc_ValueError, "Could not read FCZ file");
        return NULL;
    }
    std::vector<AtomCoordinate> atomCoordinates;
    flag = compRes.decompress(atomCoordinates);
    if (flag != 0) {
        PyErr_SetString(PyExc_ValueError, "Could not decompress FCZ file");
        return NULL;
    }

    std::vector<float3d> coordsVector = extractCoordinates(atomCoordinates);

    // Output
    PyObject* dict = getPyDictFromFoldcomp(&compRes, coordsVector);
    if (dict == NULL) {
        return NULL;
    }
    // Return dictionary
    return dict;
}

// 02. Extract information starting from PDB
PyObject* getDataFromPDB(const std::string& pdb_input) {
    std::vector<AtomCoordinate> atomCoordinates;
    // parse ATOM lines from PDB file into atomCoordinates
    std::istringstream iss(pdb_input);
    std::string line;
    // Read PDB string
    while (std::getline(iss, line)) {
        if (line.substr(0, 4) == "ATOM") {
            atomCoordinates.emplace_back(
                trim(line.substr(12, 4)), // atom
                trim(line.substr(17, 3)), // residue
                line.substr(21, 1), // chain
                std::stoi(line.substr(6, 5)), // atom_index
                std::stoi(line.substr(22, 4)), // residue_index
                std::stof(line.substr(30, 8)), std::stof(line.substr(38, 8)), std::stof(line.substr(46, 8)), // coordinates
                std::stof(line.substr(54, 6)), // occupancy
                std::stof(line.substr(60, 6)) // tempFactor
            );
        }
    }
    if (atomCoordinates.size() == 0) {
        PyErr_SetString(PyExc_ValueError, "No ATOM lines found in PDB file");
        return NULL;
    }

    // compress
    Foldcomp compRes;
    compRes.compress(atomCoordinates);

    std::vector<float3d> coordsVector = extractCoordinates(atomCoordinates);

    // Output
    PyObject* dict = getPyDictFromFoldcomp(&compRes, coordsVector);
    if (dict == NULL) {
        return NULL;
    }
    return dict;
}

static PyObject* foldcomp_get_data(PyObject* /* self */, PyObject* args, PyObject* kwargs) {
    const char* input;
    size_t input_size;
    static const char* kwlist[] = { "input", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#", (char**)kwlist, &input, &input_size)) {
        return NULL;
    }
    // Check input
    if (input_size == 0) {
        PyErr_SetString(PyExc_ValueError, "Input is empty");
        return NULL;
    }
    // Check the first 4 bytes of the input and if they are "FCMP" then it is a FCZ file
    if (input_size >= 4 && strncmp(input, "FCMP", 4) == 0) {
        return getDataFromFCZ(input, input_size);
    } else if (input_size >= 4) {
        std::string pdb_input(input, input_size);
        return getDataFromPDB(pdb_input);
    } else {
        PyErr_SetString(PyExc_ValueError, "Input is not a FCZ file or PDB file");
        return NULL;
    }
}

// Method definitions
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wcast-function-type"
static PyMethodDef foldcomp_methods[] = {
    // {"compress", foldcomp_compress, METH_VARARGS, "Compress a PDB file."},
    {"decompress", foldcomp_decompress, METH_VARARGS, "Decompress FCZ content to PDB."},
    {"compress", (PyCFunction)foldcomp_compress, METH_VARARGS | METH_KEYWORDS, "Compress PDB content to FCZ."},
    {"open", (PyCFunction)foldcomp_open, METH_VARARGS | METH_KEYWORDS, "Open a Foldcomp database."},
    {"get_data", (PyCFunction)foldcomp_get_data, METH_VARARGS | METH_KEYWORDS, "Get data from FCZ or PDB content."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};
#pragma GCC diagnostic pop
// Module definition
static struct PyModuleDef foldcomp_module_def = {
    PyModuleDef_HEAD_INIT,
    "foldcomp", /* m_name */
    NULL, /* m_doc */
    -1, /* m_size */
    foldcomp_methods, /* m_methods */
    0, /* m_slots */
    0, /* m_traverse */
    0, /* m_clear */
    0, /* m_free */
};
// Module initialization
PyMODINIT_FUNC PyInit_foldcomp(void) {
    if (PyType_Ready(&FoldcompDatabaseType) < 0) {
        return NULL;
    }

    PyObject *m = PyModule_Create(&foldcomp_module_def);
    if (m == NULL) {
        return NULL;
    }

    FoldcompError = PyErr_NewException("foldcomp.error", NULL, NULL);
    Py_XINCREF(FoldcompError);
    if (PyModule_AddObject(m, "error", FoldcompError) < 0) {
        goto clean_err;
    }

    Py_INCREF(&FoldcompDatabaseType);
    if (PyModule_AddObject(m, "FoldcompDatabase", (PyObject *)&FoldcompDatabaseType) < 0) {
        goto clean_db;
    }

    return m;

clean_db:
    Py_DECREF(&FoldcompDatabaseType);

clean_err:
    Py_XDECREF(FoldcompError);
    Py_CLEAR(FoldcompError);

    Py_DECREF(m);

    return NULL;
}
