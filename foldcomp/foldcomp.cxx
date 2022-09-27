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
#include "dbreader.h"

static PyObject *FoldcompError;

// int compress(const char* buffer, char* output, int* output_size) {
//     return 0;
// }

typedef struct {
    PyObject_HEAD
    PyObject* uniprot_ids;
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
    if (db->uniprot_ids != NULL) {
        return PySequence_Length(db->uniprot_ids);
    }
    return (Py_ssize_t)reader_get_size(db->memory_handle);
}

// FoldcompDatabase_sq_item
static PyObject* FoldcompDatabase_sq_item(PyObject* self, Py_ssize_t index) {
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;

    const char* data;
    size_t length;
    if (db->uniprot_ids != NULL) {
        if (index >= PySequence_Length(db->uniprot_ids)) {
           return NULL;
        }
        PyObject* item = PySequence_GetItem(db->uniprot_ids, index);
        // get string representation of id as c string
        uint32_t key = reader_lookup_entry(db->memory_handle, PyUnicode_AsUTF8(item));
        if (key == UINT32_MAX) {
            PyErr_SetString(PyExc_KeyError, "Could not find key in database.");
            return NULL;
        }
        int64_t id = reader_get_id(db->memory_handle, key);
        if (id == -1) {
            return NULL;
        }
        data = reader_get_data(db->memory_handle, id);
        length = std::max(reader_get_length(db->memory_handle, id), (int64_t)1) - (int64_t)1;
    } else {
        if (index >= reader_get_size(db->memory_handle)) {
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
        return Py_BuildValue("(s,O)", name.c_str(), PyBytes_FromStringAndSize(oss.str().c_str(), oss.str().size()));
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
         return NULL;
    }
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;
    Py_XDECREF(db->uniprot_ids);
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

    return Py_BuildValue("(s,O)", name.c_str(), PyBytes_FromStringAndSize(oss.str().c_str(), oss.str().size()));
}

PyTypeObject* pathType = NULL;

static PyObject *foldcomp_open(PyObject* /* self */, PyObject* args, PyObject* kwargs) {
    PyObject* path;
    PyObject* uniprot_ids = NULL;
    PyObject* decompress = NULL;
    static const char *kwlist[] = {"path", "uniprot_ids", "decompress", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|$OO", const_cast<char**>(kwlist), &path, &uniprot_ids, &decompress)) {
        return NULL;
    }

    // get path-like object's path
    PyObject* pathStr = PyObject_CallMethod(path, "as_posix", NULL);
    if (pathStr == NULL) {
        return NULL;
    }
    const char* pathCStr = PyUnicode_AsUTF8(pathStr);
    if (pathCStr == NULL) {
        return NULL;
    }

    if (uniprot_ids != NULL && !PyList_Check(uniprot_ids)) {
        PyErr_SetString(PyExc_TypeError, "uniprot_ids must be a list.");
        return NULL;
    }

    if (decompress != NULL && !PyBool_Check(decompress)) {
        PyErr_SetString(PyExc_TypeError, "decompress must be a boolean");
        return NULL;
    }

    std::string dbname(pathCStr);
    std::string index = dbname + ".index";

    FoldcompDatabaseObject *obj = PyObject_New(FoldcompDatabaseObject, &FoldcompDatabaseType);
    if (obj == NULL) {
        return NULL;
    }

    obj->uniprot_ids = NULL;
    int mode = DB_READER_USE_DATA;
    if (uniprot_ids != NULL && PySequence_Length(uniprot_ids) > 0) {
        mode |= DB_READER_USE_LOOKUP;
        obj->uniprot_ids = uniprot_ids;
        Py_INCREF(obj->uniprot_ids);
    }

    if (decompress == NULL) {
        obj->decompress = true;
    } else {
        obj->decompress = PyObject_IsTrue(decompress);
    }

    obj->memory_handle = make_reader(dbname.c_str(), index.c_str(), mode);
    return (PyObject*)obj;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wcast-function-type"
static PyMethodDef foldcomp_methods[] = {
    // {"compress", foldcomp_compress, METH_VARARGS, "Compress a PDB file."},
    {"decompress", foldcomp_decompress, METH_VARARGS, "Decompress a PDB file."},
    {"open", (PyCFunction)foldcomp_open, METH_VARARGS | METH_KEYWORDS, "Open a Foldcomp database."},

    {NULL, NULL, 0, NULL} /* Sentinel */
};
#pragma GCC diagnostic pop

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
