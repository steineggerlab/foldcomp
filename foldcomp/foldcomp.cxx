#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>
#include <sstream> // IWYU pragma: keep
#include <unordered_map>

#include "atom_coordinate.h"
#include "foldcomp.h"
#include "database_reader.h"
#include "utility.h"

static PyObject *FoldcompError;

enum class OutputFormat {
    PDB,
    MMCIF
};

typedef struct {
    PyObject_HEAD
    std::vector<int64_t>* user_indices;
    std::vector<std::string>* merged_names;
    std::vector<std::vector<int64_t>>* merged_source_indices;
    bool merge_fragments;
    bool decompress;
    OutputFormat output_format;
    void* memory_handle;
} FoldcompDatabaseObject;

int decompressToAtoms(const char* input, size_t input_size, bool use_alt_order, std::vector<AtomCoordinate>& atomCoordinates, std::string& name);
int decompress(const char* input, size_t input_size, bool use_alt_order, OutputFormat format, std::ostream& oss, std::string& name);
int parsePDB(const std::string& pdb_input, std::vector<AtomCoordinate>& atomCoordinates);
static PyObject* FoldcompDatabase_close(PyObject* self);
static PyObject* FoldcompDatabase_enter(PyObject* self);
static PyObject* FoldcompDatabase_exit(PyObject* self, PyObject* args);
static PyObject* FoldcompDatabase_source_indices(PyObject* self, PyObject* args);
void FoldcompDatabase_release_resources(FoldcompDatabaseObject* db);
PyObject* vectorToList_Int64(const std::vector<int64_t>& data);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wcast-function-type"
static PyMethodDef FoldcompDatabase_methods[] = {
    {"close", (PyCFunction)FoldcompDatabase_close, METH_NOARGS, "Close the database."},
    {"source_indices", (PyCFunction)FoldcompDatabase_source_indices, METH_VARARGS, "Get source fragment indices for a merged entry."},
    {"__enter__", (PyCFunction)FoldcompDatabase_enter, METH_NOARGS, "Enter the runtime context related to this object."},
    {"__exit__", (PyCFunction)FoldcompDatabase_exit, METH_VARARGS, "Exit the runtime context related to this object."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};
#pragma GCC diagnostic pop

// FoldcompDatabase_sq_length
static Py_ssize_t FoldcompDatabase_sq_length(PyObject* self) {
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;
    if (db->merge_fragments && db->merged_names != NULL) {
        return (Py_ssize_t)db->merged_names->size();
    }
    if (db->user_indices != NULL) {
        return db->user_indices->size();
    }
    return (Py_ssize_t)reader_get_size(db->memory_handle);
}

static void appendPDBFragment(std::ostringstream& out, const std::string& pdb, bool include_title) {
    std::istringstream iss(pdb);
    std::string line;
    while (std::getline(iss, line)) {
        if (!include_title && stringStartsWith("TITLE", line)) {
            continue;
        }
        out << line << "\n";
    }
}

static std::string toLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return (char)std::tolower(c);
    });
    return value;
}

static bool parseOutputFormat(const char* format, OutputFormat& output_format) {
    if (format == NULL) {
        output_format = OutputFormat::PDB;
        return true;
    }
    std::string format_lower = toLower(std::string(format));
    if (format_lower == "pdb") {
        output_format = OutputFormat::PDB;
        return true;
    }
    if (format_lower == "mmcif" || format_lower == "cif") {
        output_format = OutputFormat::MMCIF;
        return true;
    }
    return false;
}

static void writeAtomsByFormat(std::vector<AtomCoordinate>& atomCoordinates, const std::string& name, OutputFormat output_format, std::ostream& oss) {
    if (output_format == OutputFormat::MMCIF) {
        writeAtomCoordinatesToMMCIF(atomCoordinates, name, oss);
    } else {
        writeAtomCoordinatesToPDB(atomCoordinates, name, oss);
    }
}

// FoldcompDatabase_sq_item
static PyObject* FoldcompDatabase_sq_item(PyObject* self, Py_ssize_t index) {
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;

    if (db->merge_fragments) {
        if (!db->decompress) {
            PyErr_SetString(PyExc_TypeError, "merge_fragments requires decompress=True");
            return NULL;
        }
        if (db->merged_names == NULL || db->merged_source_indices == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "merged database state is not initialized");
            return NULL;
        }
        if (index < 0 || index >= (Py_ssize_t)db->merged_names->size()) {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }

        const std::string& merged_name = db->merged_names->at(index);
        const std::vector<int64_t>& source_ids = db->merged_source_indices->at(index);
        if (db->output_format == OutputFormat::PDB) {
            std::ostringstream merged_oss;
            bool include_title = true;
            for (int64_t source_id : source_ids) {
                const char* data = reader_get_data(db->memory_handle, source_id);
                size_t length = std::max(reader_get_length(db->memory_handle, source_id), (int64_t)1) - (int64_t)1;
                if (data == NULL) {
                    PyErr_SetString(FoldcompError, "Failed to read source fragment from database.");
                    return NULL;
                }
                std::ostringstream fragment_oss;
                std::string fragment_name;
                int err = decompress(data, length, false, OutputFormat::PDB, fragment_oss, fragment_name);
                if (err != 0) {
                    std::string err_msg = "Error decompressing: " + fragment_name;
                    PyErr_SetString(FoldcompError, err_msg.c_str());
                    return NULL;
                }
                std::string fragment_pdb = fragment_oss.str();
                appendPDBFragment(merged_oss, fragment_pdb, include_title);
                include_title = false;
            }
            std::string merged_pdb = merged_oss.str();
            PyObject* pdb = PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, merged_pdb.c_str(), merged_pdb.size());
            PyObject* result = Py_BuildValue("(s,O)", merged_name.c_str(), pdb);
            Py_DECREF(pdb);
            return result;
        }

        std::vector<AtomCoordinate> merged_atoms;
        for (int64_t source_id : source_ids) {
            const char* data = reader_get_data(db->memory_handle, source_id);
            size_t length = std::max(reader_get_length(db->memory_handle, source_id), (int64_t)1) - (int64_t)1;
            if (data == NULL) {
                PyErr_SetString(FoldcompError, "Failed to read source fragment from database.");
                return NULL;
            }
            std::vector<AtomCoordinate> fragment_atoms;
            std::string fragment_name;
            int err = decompressToAtoms(data, length, false, fragment_atoms, fragment_name);
            if (err != 0) {
                std::string err_msg = "Error decompressing: " + fragment_name;
                PyErr_SetString(FoldcompError, err_msg.c_str());
                return NULL;
            }
            merged_atoms.insert(merged_atoms.end(), fragment_atoms.begin(), fragment_atoms.end());
        }

        std::ostringstream merged_oss;
        writeAtomsByFormat(merged_atoms, merged_name, db->output_format, merged_oss);
        std::string merged_text = merged_oss.str();
        PyObject* structure = PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, merged_text.c_str(), merged_text.size());
        PyObject* result = Py_BuildValue("(s,O)", merged_name.c_str(), structure);
        Py_DECREF(structure);
        return result;
    }

    const char* data = NULL;
    size_t length = 0;
    int64_t id = -1;
    if (db->user_indices != NULL) {
        if (index >= (Py_ssize_t)db->user_indices->size()) {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }
        id = db->user_indices->at(index);
        data = reader_get_data(db->memory_handle, id);
        length = std::max(reader_get_length(db->memory_handle, id), (int64_t)1) - (int64_t)1;
    } else {
        if (index >= (Py_ssize_t)reader_get_size(db->memory_handle)) {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }
        data = reader_get_data(db->memory_handle, index);
        length = std::max(reader_get_length(db->memory_handle, index), (int64_t)1) - (int64_t)1;
    }
    if (db->decompress) {
        std::ostringstream oss;
        std::string name;
        int err = decompress(data, length, false, db->output_format, oss, name);
        if (err != 0) {
            std::string err_msg = "Error decompressing: " + name;
            PyErr_SetString(FoldcompError, err_msg.c_str());
            return NULL;
        }
        std::string structure_text = oss.str();
        PyObject* pdb = PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, structure_text.c_str(), structure_text.size());
        PyObject* result = Py_BuildValue("(s,O)", name.c_str(), pdb);
        Py_DECREF(pdb);
        return result;
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
    FoldcompDatabase_release_resources(db);
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

static PyObject* FoldcompDatabase_source_indices(PyObject* self, PyObject* args) {
    if (!PyObject_TypeCheck(self, &FoldcompDatabaseType)) {
        PyErr_SetString(PyExc_TypeError, "Expected FoldcompDatabase object.");
        return NULL;
    }
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;
    if (!db->merge_fragments || db->merged_source_indices == NULL) {
        PyErr_SetString(PyExc_TypeError, "source_indices is available only when merge_fragments=True.");
        return NULL;
    }
    Py_ssize_t index;
    if (!PyArg_ParseTuple(args, "n", &index)) {
        return NULL;
    }
    if (index < 0 || index >= (Py_ssize_t)db->merged_source_indices->size()) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return vectorToList_Int64(db->merged_source_indices->at(index));
}

void FoldcompDatabase_release_resources(FoldcompDatabaseObject* db) {
    if (db == NULL) {
        return;
    }
    if (db->memory_handle != NULL) {
        free_reader(db->memory_handle);
        db->memory_handle = NULL;
    }
    if (db->user_indices != NULL) {
        delete db->user_indices;
        db->user_indices = NULL;
    }
    if (db->merged_names != NULL) {
        delete db->merged_names;
        db->merged_names = NULL;
    }
    if (db->merged_source_indices != NULL) {
        delete db->merged_source_indices;
        db->merged_source_indices = NULL;
    }
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
int decompressToAtoms(const char* input, size_t input_size, bool use_alt_order, std::vector<AtomCoordinate>& atomCoordinates, std::string& name) {
    OneShotReadBuf buf((char*)input, input_size);
    std::istream istr(&buf);

    std::ios_base::iostate cout_state = std::cout.rdstate();
    std::cout.setstate(std::ios_base::failbit);
    Foldcomp compRes;
    int flag = compRes.read(istr);
    if (flag != 0) {
        std::cout.clear(cout_state);
        return 1;
    }
    compRes.useAltAtomOrder = use_alt_order;
    flag = compRes.decompress(atomCoordinates);
    if (flag != 0) {
        std::cout.clear(cout_state);
        return 1;
    }
    std::cout.clear(cout_state);

    name = compRes.strTitle;

    return 0;
}

int decompress(const char* input, size_t input_size, bool use_alt_order, OutputFormat format, std::ostream& oss, std::string& name) {
    std::vector<AtomCoordinate> atomCoordinates;
    int flag = decompressToAtoms(input, input_size, use_alt_order, atomCoordinates, name);
    if (flag != 0) {
        return flag;
    }
    writeAtomsByFormat(atomCoordinates, name, format, oss);
    return 0;
}
// Python binding for decompress
static PyObject *foldcomp_decompress(PyObject* /* self */, PyObject *args, PyObject* kwargs) {
    // Unpack a string from the arguments
    const char *strArg;
    Py_ssize_t strSize;
    const char* format = "pdb";
    static const char *kwlist[] = {"input", "format", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "y#|$s", const_cast<char**>(kwlist), &strArg, &strSize, &format)) {
        return NULL;
    }
    OutputFormat output_format;
    if (!parseOutputFormat(format, output_format)) {
        PyErr_SetString(PyExc_ValueError, "format must be one of: 'pdb', 'mmcif', 'cif'");
        return NULL;
    }

    std::ostringstream oss;
    std::string name;
    int err = decompress(strArg, strSize, false, output_format, oss, name);
    if (err != 0) {
        PyErr_SetString(FoldcompError, "Error decompressing.");
        return NULL;
    }

    std::string structure_text = oss.str();
    return Py_BuildValue("(s,O)", name.c_str(), PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, structure_text.c_str(), structure_text.size()));
}

std::string trim(const std::string& str, const std::string& whitespace = " \t") {
    const std::string::size_type strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const std::string::size_type strEnd = str.find_last_not_of(whitespace);
    const std::string::size_type strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

int parsePDB(const std::string& pdb_input, std::vector<AtomCoordinate>& atomCoordinates) {
    atomCoordinates.clear();
    std::istringstream iss(pdb_input);
    std::string line;
    std::string chain = "";
    while (std::getline(iss, line)) {
        if (line.substr(0, 4) == "ATOM") {
            if (chain == "") {
                chain = line.substr(21, 1);
            }
            if (line.substr(21, 1) != chain) {
                return 2; // FLAG 2: multiple chains
            }
            atomCoordinates.emplace_back(
                trim(line.substr(12, 4)), // atom
                trim(line.substr(17, 3)), // residue
                chain, // chain
                std::stoi(line.substr(6,  5)), // atom_index
                std::stoi(line.substr(22, 4)), // residue_index
                std::stof(line.substr(30, 8)), std::stof(line.substr(38, 8)), std::stof(line.substr(46, 8)), // coordinates
                std::stof(line.substr(54, 6)), // occupancy
                std::stof(line.substr(60, 6)) // tempFactor
            );
        }
    }
    if (atomCoordinates.size() == 0) {
        return 1; // FLAG 1: no ATOM lines
    }
    return 0;
}

// Compress
int compress(const std::string& name, const std::string& pdb_input, std::ostream& oss, int anchor_residue_threshold) {
    std::vector<AtomCoordinate> atomCoordinates;
    int parse_flag = parsePDB(pdb_input, atomCoordinates);
    if (parse_flag != 0) {
        return parse_flag;
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
    PyObject* split = NULL;
    static const char *kwlist[] = {"name", "pdb_content", "anchor_residue_threshold", "split", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|$OO", const_cast<char**>(kwlist), &name, &pdb_input, &anchor_residue_threshold, &split)) {
        return NULL;
    }

    if (anchor_residue_threshold != NULL && !PyLong_Check(anchor_residue_threshold)) {
        PyErr_SetString(PyExc_TypeError, "anchor_residue_threshold must be an integer");
        return NULL;
    }
    if (split != NULL && !PyBool_Check(split)) {
        PyErr_SetString(PyExc_TypeError, "split must be a boolean");
        return NULL;
    }

    int threshold = DEFAULT_ANCHOR_THRESHOLD;
    if (anchor_residue_threshold != NULL) {
        threshold = PyLong_AsLong(anchor_residue_threshold);
    }
    bool split_flag = split != NULL && PyObject_IsTrue(split);

    if (split_flag) {
        std::string pdb_content(pdb_input);
        std::vector<std::pair<std::string, std::string>> chains;
        std::istringstream chain_stream(pdb_content);
        std::string line;
        std::string current_chain = "";
        std::string current_chain_chunk = "";
        while (std::getline(chain_stream, line)) {
            if (!stringStartsWith("ATOM", line)) {
                continue;
            }
            std::string chain = line.substr(21, 1);
            if (current_chain.empty()) {
                current_chain = chain;
            } else if (chain != current_chain) {
                if (!current_chain_chunk.empty()) {
                    chains.emplace_back(current_chain, current_chain_chunk);
                }
                current_chain = chain;
                current_chain_chunk.clear();
            }
            current_chain_chunk += line + "\n";
        }
        if (!current_chain_chunk.empty()) {
            chains.emplace_back(current_chain, current_chain_chunk);
        }
        if (chains.size() == 0) {
            PyErr_SetString(FoldcompError, "No ATOM lines found");
            return NULL;
        }

        std::string base = baseName(name);
        std::pair<std::string, std::string> output_parts = getFileParts(base);
        std::string output_base = output_parts.first;

        PyObject* outputs = PyList_New(0);
        if (outputs == NULL) {
            return NULL;
        }

        for (size_t i = 0; i < chains.size(); i++) {
            std::vector<std::string> fragments;
            std::istringstream fragment_stream(chains[i].second);
            std::string fragment;
            int prev_n_res_idx = 0;
            bool has_prev_n = false;
            while (std::getline(fragment_stream, line)) {
                if (!stringStartsWith("ATOM", line)) {
                    continue;
                }
                std::string atom = trim(line.substr(12, 4));
                if (atom == "N") {
                    int curr_res_idx = 0;
                    bool parsed = true;
                    try {
                        curr_res_idx = std::stoi(line.substr(22, 4));
                    } catch (...) {
                        parsed = false;
                    }
                    if (parsed && has_prev_n && curr_res_idx - prev_n_res_idx > 1) {
                        if (!fragment.empty()) {
                            fragments.push_back(fragment);
                            fragment.clear();
                        }
                    }
                    if (parsed) {
                        prev_n_res_idx = curr_res_idx;
                        has_prev_n = true;
                    }
                }
                fragment += line + "\n";
            }
            if (!fragment.empty()) {
                fragments.push_back(fragment);
            }

            for (size_t j = 0; j < fragments.size(); j++) {
                std::ostringstream oss;
                int flag = compress(name, fragments[j], oss, threshold);
                if (flag != 0) {
                    continue;
                }

                std::string chunk_name = output_base;
                if (chains.size() > 1) {
                    chunk_name += chains[i].first;
                }
                if (fragments.size() > 1) {
                    chunk_name += "_" + std::to_string(j);
                }
                chunk_name += ".fcz";

                std::string payload = oss.str();
                PyObject* py_payload = PyBytes_FromStringAndSize(payload.c_str(), payload.size());
                if (py_payload == NULL) {
                    Py_DECREF(outputs);
                    return NULL;
                }
                PyObject* py_item = Py_BuildValue("(s,O)", chunk_name.c_str(), py_payload);
                Py_DECREF(py_payload);
                if (py_item == NULL) {
                    Py_DECREF(outputs);
                    return NULL;
                }
                if (PyList_Append(outputs, py_item) != 0) {
                    Py_DECREF(py_item);
                    Py_DECREF(outputs);
                    return NULL;
                }
                Py_DECREF(py_item);
            }
        }
        return outputs;
    }

    std::ostringstream oss;
    int flag = compress(name, pdb_input, oss, threshold);
    if (flag == 1) {
        PyErr_SetString(FoldcompError, "No ATOM lines found");
        return NULL;
    } else if (flag == 2) {
        PyErr_SetString(FoldcompError, "Multiple chains found. Please provide a single chain using 'foldcomp.split_pdb_by_chain'");
        return NULL;
    } else if (flag != 0) {
        PyErr_SetString(FoldcompError, "Error compressing");
        return NULL;
    }

    std::string compressed_bytes = oss.str();
    return PyBytes_FromStringAndSize(compressed_bytes.c_str(), compressed_bytes.length());
}


PyTypeObject* pathType = NULL;

static PyObject *foldcomp_open(PyObject* /* self */, PyObject* args, PyObject* kwargs) {
    PyObject* path;
    PyObject* user_ids = NULL;
    PyObject* decompress = NULL;
    PyObject* err_on_missing = NULL; // Raise an error if the file is missing. Default: False
    PyObject* merge_fragments = NULL;
    const char* format = "pdb";

    static const char *kwlist[] = {"path", "ids", "decompress", "err_on_missing", "merge_fragments", "format", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&|$OOOOs", const_cast<char**>(kwlist), PyUnicode_FSConverter, &path, &user_ids, &decompress, &err_on_missing, &merge_fragments, &format)) {
        return NULL;
    }
    if (path == NULL) {
        PyErr_SetString(PyExc_TypeError, "path must be a path-like object");
        return NULL;
    }
    const char* pathCStr = PyBytes_AS_STRING(path);
    if (pathCStr == NULL) {
        Py_DECREF(path);
        PyErr_SetString(PyExc_TypeError, "path must be a path-like object");
        return NULL;
    }

    if (user_ids != NULL && !PyList_Check(user_ids)) {
        Py_DECREF(path);
        PyErr_SetString(PyExc_TypeError, "user_ids must be a list.");
        return NULL;
    }

    if (decompress != NULL && !PyBool_Check(decompress)) {
        Py_DECREF(path);
        PyErr_SetString(PyExc_TypeError, "decompress must be a boolean");
        return NULL;
    }

    if (err_on_missing != NULL && !PyBool_Check(err_on_missing)) {
        Py_DECREF(path);
        PyErr_SetString(PyExc_TypeError, "err_on_missing must be a boolean");
        return NULL;
    }
    if (merge_fragments != NULL && !PyBool_Check(merge_fragments)) {
        Py_DECREF(path);
        PyErr_SetString(PyExc_TypeError, "merge_fragments must be a boolean");
        return NULL;
    }
    OutputFormat output_format;
    if (!parseOutputFormat(format, output_format)) {
        Py_DECREF(path);
        PyErr_SetString(PyExc_ValueError, "format must be one of: 'pdb', 'mmcif', 'cif'");
        return NULL;
    }

    std::string dbname(pathCStr);
    std::string index = dbname + ".index";
    bool err_on_missing_flag = false;
    Py_DECREF(path);

    FoldcompDatabaseObject *obj = PyObject_New(FoldcompDatabaseObject, &FoldcompDatabaseType);
    if (obj == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for FoldcompDatabaseObject");
        return NULL;
    }
    obj->memory_handle = NULL;
    obj->user_indices = NULL;
    obj->merged_names = NULL;
    obj->merged_source_indices = NULL;
    obj->merge_fragments = false;
    obj->decompress = true;
    obj->output_format = output_format;

    int mode = DB_READER_USE_DATA;
    bool merge_fragments_flag = merge_fragments != NULL && PyObject_IsTrue(merge_fragments);
    if (merge_fragments_flag) {
        mode |= DB_READER_USE_LOOKUP_REVERSE;
    } else if (user_ids != NULL && PySequence_Length(user_ids) > 0) {
        mode |= DB_READER_USE_LOOKUP;
    }

    if (decompress == NULL) {
        obj->decompress = true;
    } else {
        obj->decompress = PyObject_IsTrue(decompress);
    }

    if (err_on_missing == NULL) {
        err_on_missing_flag = false;
    } else {
        err_on_missing_flag = PyObject_IsTrue(err_on_missing);
    }
    obj->merge_fragments = merge_fragments_flag;
    if (obj->merge_fragments && !obj->decompress) {
        FoldcompDatabase_release_resources(obj);
        Py_DECREF(obj);
        PyErr_SetString(PyExc_TypeError, "merge_fragments requires decompress=True");
        return NULL;
    }

    obj->memory_handle = make_reader(dbname.c_str(), index.c_str(), mode);
    if (obj->memory_handle == NULL) {
        FoldcompDatabase_release_resources(obj);
        Py_DECREF(obj);
        PyErr_SetString(PyExc_RuntimeError, "Failed to open Foldcomp database.");
        return NULL;
    }

    if (obj->merge_fragments) {
        obj->merged_names = new std::vector<std::string>();
        obj->merged_source_indices = new std::vector<std::vector<int64_t>>();
        std::unordered_map<std::string, size_t> group_idx_by_name;
        int64_t size = reader_get_size(obj->memory_handle);
        for (int64_t id = 0; id < size; id++) {
            uint32_t key = reader_get_key(obj->memory_handle, id);
            const char* lookup_name = reader_lookup_name_alloc(obj->memory_handle, key);
            std::string group_name;
            if (lookup_name != NULL && lookup_name[0] != '\0') {
                group_name = lookup_name;
                free((void*)lookup_name);
            } else {
                group_name = std::to_string(key);
            }

            auto it = group_idx_by_name.find(group_name);
            size_t group_idx;
            if (it == group_idx_by_name.end()) {
                group_idx = obj->merged_names->size();
                group_idx_by_name[group_name] = group_idx;
                obj->merged_names->push_back(group_name);
                obj->merged_source_indices->emplace_back();
            } else {
                group_idx = it->second;
            }
            obj->merged_source_indices->at(group_idx).push_back(id);
        }

        if (user_ids != NULL && PySequence_Length(user_ids) > 0) {
            std::vector<std::string> filtered_names;
            std::vector<std::vector<int64_t>> filtered_sources;
            size_t id_count = (size_t)PySequence_Length(user_ids);
            for (Py_ssize_t i = 0; i < (Py_ssize_t)id_count; i++) {
                PyObject* item = PySequence_GetItem(user_ids, i);
                const char* data = PyUnicode_AsUTF8(item);
                Py_DECREF(item);
                auto it = group_idx_by_name.find(data);
                if (it == group_idx_by_name.end()) {
                    std::string err_msg = "Skipping entry ";
                    err_msg += data;
                    err_msg += " which is not in the database.";
                    if (err_on_missing_flag) {
                        FoldcompDatabase_release_resources(obj);
                        Py_DECREF(obj);
                        PyErr_SetString(PyExc_KeyError, err_msg.c_str());
                        return NULL;
                    } else {
                        std::cerr << err_msg << std::endl;
                        continue;
                    }
                }
                filtered_names.push_back(data);
                filtered_sources.push_back(obj->merged_source_indices->at(it->second));
            }
            *(obj->merged_names) = filtered_names;
            *(obj->merged_source_indices) = filtered_sources;
        }
    } else if (user_ids != NULL && PySequence_Length(user_ids) > 0) {
        size_t id_count = (size_t)PySequence_Length(user_ids);
        // Reserve memory for the user indices
        obj->user_indices = new std::vector<int64_t>();
        obj->user_indices->reserve(id_count);
        for (Py_ssize_t i = 0; i < (Py_ssize_t)id_count; i++) {
            // Iterate over all entries in the database and store ids in a vector of int64_t
            PyObject* item = PySequence_GetItem(user_ids, i);
            const char* data = PyUnicode_AsUTF8(item);
            Py_DECREF(item);
            uint32_t key = reader_lookup_entry(obj->memory_handle, data);
            int64_t id = reader_get_id(obj->memory_handle, key);
            if (id == -1 || key == UINT32_MAX) {
                // Not found --> no error just
                std::string err_msg = "Skipping entry ";
                err_msg += data;
                err_msg += " which is not in the database.";
                if (err_on_missing_flag) {
                    FoldcompDatabase_release_resources(obj);
                    Py_DECREF(obj);
                    PyErr_SetString(PyExc_KeyError, err_msg.c_str());
                    return NULL;
                } else {
                    std::cerr << err_msg << std::endl;
                    continue;
                }
            }
            obj->user_indices->push_back(id);
        }
    }

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

PyObject* vectorToList_Int64(const std::vector<int64_t>& data) {
    PyObject* listObj = PyList_New(data.size());
    if (!listObj) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for list");
        return NULL;
    }
    for (size_t i = 0; i < data.size(); i++) {
        // data[i] is a int64_t, but PyLong_FromLongLong expects a long long
        // so we need to cast it without error
        PyObject* num = PyLong_FromLongLong((long long)data[i]);
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

    // Free memory
    Py_XDECREF(phi);
    Py_XDECREF(psi);
    Py_XDECREF(omega);
    Py_XDECREF(torsion_angles);
    Py_XDECREF(bond_angles);
    Py_XDECREF(residues);
    Py_XDECREF(b_factors);
    Py_XDECREF(coordinates);

    return dict;
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
    // Free memory
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
    {"decompress", (PyCFunction)foldcomp_decompress, METH_VARARGS | METH_KEYWORDS, "Decompress FCZ content to PDB or mmCIF."},
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
