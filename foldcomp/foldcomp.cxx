#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <atomic>
#include <cstdint>
#include <cstddef>
#include <exception>
#include <limits>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>
#include <utility>

#include "amino_acid.h"
#include "atom_coordinate.h"
#include "foldcomp.h"
#include "database_reader.h"
#include "structure_codec.h"
#ifdef FOLDCOMP_WITH_CUDA
#include <cuda_runtime.h>
#include <nvtx3/nvtx3.hpp>  // FOLDCOMP_CPU_NVTX + explicit ranges for host tail
#include "gpu_decompression_pipeline.h"
#include "gpu_pdb_writer.h"
#include "input_processor.h"
#endif

static PyObject *FoldcompError;

typedef struct {
    PyObject_HEAD
    std::vector<int64_t>* user_indices;
    bool decompress;
    void* memory_handle;
} FoldcompDatabaseObject;

int decompress(const char* input, size_t input_size, bool use_alt_order, std::string& output, std::string& name);
int decompress(const char* input, size_t input_size, bool use_alt_order, const std::string& format, std::string& output, std::string& name);
PyObject* getDataFromStructureText(const std::string& input, const char* format);
static PyObject* FoldcompDatabase_close(PyObject* self);
static PyObject* FoldcompDatabase_enter(PyObject* self);
static PyObject* FoldcompDatabase_exit(PyObject* self, PyObject* args);
static void FoldcompDatabase_dealloc(PyObject* self);
PyObject* vectorToList_Int64(const std::vector<int64_t>& data);

namespace {

bool regionNeedsRawFallback(const tcb::span<AtomCoordinate>& atoms);
bool shouldStoreSmallMixedFragmentAsRaw(
    const tcb::span<AtomCoordinate>& chainSpan,
    const std::vector<BackboneRegion>& regions,
    const std::vector<bool>& regionNeedsRaw
);

static int encodeStructureToFoldcompContainer(
    const std::string& title,
    const char* data,
    size_t size,
    const char* format,
    int anchorResidueThreshold,
    float maxBackboneRmsd,
    std::string& output
) {
    std::vector<AtomCoordinate> atomCoordinates;
    int status = PARSE_PDB_OK;
    if (!parseStructureAtoms(data, size, false, atomCoordinates, status, nullptr, format)) {
        return status;
    }

    std::vector<ContainerFragment> encodedFragments;
    std::vector<std::pair<size_t, size_t>> chainIndices = identifyChains(atomCoordinates);
    for (const auto& chainRegion : chainIndices) {
        std::vector<std::pair<size_t, size_t>> fragmentIndices = identifyDiscontinousResInd(
            atomCoordinates, chainRegion.first, chainRegion.second
        );
        if (fragmentIndices.empty()) {
            fragmentIndices.push_back(chainRegion);
        }
        for (const auto& fragment : fragmentIndices) {
            tcb::span<AtomCoordinate> chainSpan(
                atomCoordinates.data() + fragment.first,
                fragment.second - fragment.first
            );
            if (chainSpan.empty()) {
                continue;
            }
            std::vector<BackboneRegion> regions = identifyBackboneRegions(chainSpan);
            if (regions.empty()) {
                continue;
            }
            std::vector<bool> regionNeedsRaw(regions.size(), false);
            for (size_t regionIndex = 0; regionIndex < regions.size(); regionIndex++) {
                if (!regions[regionIndex].encodable) {
                    continue;
                }
                tcb::span<AtomCoordinate> regionSpan(
                    chainSpan.data() + regions[regionIndex].start,
                    regions[regionIndex].end - regions[regionIndex].start
                );
                regionNeedsRaw[regionIndex] = regionNeedsRawFallback(regionSpan);
            }
            if (shouldStoreSmallMixedFragmentAsRaw(chainSpan, regions, regionNeedsRaw)) {
                ContainerFragment containerFragment;
                containerFragment.kind = CONTAINER_FRAGMENT_KIND_RAW_ATOMS;
                containerFragment.model = chainSpan.front().model;
                containerFragment.chain = chainSpan.front().chain;
                if (!serializeAtomCoordinates(
                        tcb::span<const AtomCoordinate>(chainSpan.data(), chainSpan.size()),
                        containerFragment.payload)) {
                    return PARSE_PDB_INVALID_FORMAT;
                }
                encodedFragments.push_back(std::move(containerFragment));
                continue;
            }
            for (size_t regionIndex = 0; regionIndex < regions.size(); regionIndex++) {
                const auto& region = regions[regionIndex];
                tcb::span<AtomCoordinate> regionSpan(
                    chainSpan.data() + region.start,
                    region.end - region.start
                );
                ContainerFragment containerFragment;
                containerFragment.model = chainSpan[region.start].model;
                containerFragment.chain = chainSpan[region.start].chain;
                bool useFoldcompEncoding = region.encodable && !regionNeedsRaw[regionIndex];
                if (useFoldcompEncoding) {
                    Foldcomp compRes;
                    compRes.strTitle = title;
                    compRes.anchorThreshold = anchorResidueThreshold;
                    std::vector<BackboneChain> compData = compRes.compress(regionSpan);
                    if (compData.empty()) {
                        continue;
                    }
                    containerFragment.kind = CONTAINER_FRAGMENT_KIND_FCZ;
                    if (compRes.writeString(containerFragment.payload) != 0) {
                        return PARSE_PDB_INVALID_FORMAT;
                    }
                    if (compRes.exceedsBackboneRmsdThreshold(maxBackboneRmsd)) {
                        containerFragment.kind = CONTAINER_FRAGMENT_KIND_RAW_ATOMS;
                        containerFragment.payload.clear();
                        if (!serializeAtomCoordinates(
                                tcb::span<const AtomCoordinate>(regionSpan.data(), regionSpan.size()),
                                containerFragment.payload)) {
                            return PARSE_PDB_INVALID_FORMAT;
                        }
                    }
                } else {
                    containerFragment.kind = CONTAINER_FRAGMENT_KIND_RAW_ATOMS;
                    if (!serializeAtomCoordinates(
                            tcb::span<const AtomCoordinate>(regionSpan.data(), regionSpan.size()),
                            containerFragment.payload)) {
                        return PARSE_PDB_INVALID_FORMAT;
                    }
                }
                encodedFragments.push_back(std::move(containerFragment));
            }
        }
    }

    if (encodedFragments.empty()) {
        return PARSE_PDB_NO_ATOM;
    }

    bool useContainer = encodedFragments.size() > 1 ||
                        encodedFragments[0].kind != CONTAINER_FRAGMENT_KIND_FCZ ||
                        encodedFragments[0].chain.size() != 1 ||
                        encodedFragments[0].model != 1;
    if (useContainer) {
        if (!writeContainerToString(output, title, encodedFragments)) {
            return PARSE_PDB_INVALID_FORMAT;
        }
    } else {
        output = encodedFragments[0].payload;
    }
    return PARSE_PDB_OK;
}

struct CoutStateGuard {
    std::ios::iostate state;
    CoutStateGuard(): state(std::cout.rdstate()) {}
    ~CoutStateGuard() {
        std::cout.clear(state);
    }
};

bool residueNeedsRawFallback(const tcb::span<const AtomCoordinate>& residueAtoms) {
    if (residueAtoms.empty()) {
        return false;
    }
    auto aaIt = Foldcomp::AAS.find(residueAtoms[0].residue);
    if (aaIt == Foldcomp::AAS.end()) {
        return true;
    }

    int matchedCanonicalAtoms = 0;
    bool hasOxt = false;
    for (const auto& atom : residueAtoms) {
        if (atom.altloc != ' ' && atom.altloc != '\0') {
            return true;
        }
        if (atom.insertion_code != ' ' && atom.insertion_code != '\0') {
            return true;
        }
        if (atom.atom == "OXT") {
            if (hasOxt) {
                return true;
            }
            hasOxt = true;
            continue;
        }
        bool found = false;
        for (const auto& canonicalAtom : aaIt->second.atoms) {
            if (atom.atom == canonicalAtom) {
                matchedCanonicalAtoms++;
                found = true;
                break;
            }
        }
        if (!found) {
            return true;
        }
    }
    int expectedAtomCount = static_cast<int>(aaIt->second.atoms.size()) + (hasOxt ? 1 : 0);
    return static_cast<int>(residueAtoms.size()) != expectedAtomCount ||
           matchedCanonicalAtoms != static_cast<int>(aaIt->second.atoms.size());
}

bool regionNeedsRawFallback(const tcb::span<AtomCoordinate>& atoms) {
    size_t residuesWithOxt = 0;
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) || startsNewResidue(atoms[i], atoms[i - 1]);
        if (!endOfResidue) {
            continue;
        }
        bool hasOxt = false;
        for (size_t j = residueStart; j < i; j++) {
            const auto& atom = atoms[j];
            if (atom.atom == "OXT") {
                hasOxt = true;
                break;
            }
        }
        if (hasOxt) {
            residuesWithOxt++;
            if (residuesWithOxt > 1 || i != atoms.size()) {
                return true;
            }
        }
        if (residueNeedsRawFallback(
                tcb::span<const AtomCoordinate>(atoms.data() + residueStart, i - residueStart))) {
            return true;
        }
        residueStart = i;
    }
    return false;
}

size_t countResidues(const tcb::span<AtomCoordinate>& atoms) {
    if (atoms.empty()) {
        return 0;
    }
    size_t residueCount = 1;
    for (size_t i = 1; i < atoms.size(); i++) {
        if (startsNewResidue(atoms[i], atoms[i - 1])) {
            residueCount++;
        }
    }
    return residueCount;
}

bool shouldStoreSmallMixedFragmentAsRaw(
    const tcb::span<AtomCoordinate>& chainSpan,
    const std::vector<BackboneRegion>& regions,
    const std::vector<bool>& regionNeedsRaw
) {
    if (regions.size() <= 1) {
        return false;
    }
    if (countResidues(chainSpan) > 24) {
        return false;
    }
    for (size_t i = 0; i < regions.size(); i++) {
        if (!regions[i].encodable || regionNeedsRaw[i]) {
            return true;
        }
    }
    return false;
}

void releaseFoldcompDatabase(FoldcompDatabaseObject* db) {
    if (db->memory_handle != NULL) {
        free_reader(db->memory_handle);
        db->memory_handle = NULL;
    }
    delete db->user_indices;
    db->user_indices = NULL;
}

}

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
    if (db->memory_handle == NULL) {
        PyErr_SetString(PyExc_ValueError, "database is closed");
        return -1;
    }
    if (db->user_indices != NULL) {
        return db->user_indices->size();
    }
    return (Py_ssize_t)reader_get_size(db->memory_handle);
}

// FoldcompDatabase_sq_item
static PyObject* FoldcompDatabase_sq_item(PyObject* self, Py_ssize_t index) {
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;
    if (db->memory_handle == NULL) {
        PyErr_SetString(PyExc_ValueError, "database is closed");
        return NULL;
    }
    if (index < 0) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }

    const char* data;
    size_t length;
    int64_t id;
    if (db->user_indices != NULL) {
        if (index >= (Py_ssize_t)db->user_indices->size()) {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }
        id = db->user_indices->at(index);
        data = reader_get_data(db->memory_handle, id);
        length = reader_get_length(db->memory_handle, id);
    } else {
        if (index >= (Py_ssize_t)reader_get_size(db->memory_handle)) {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }
        data = reader_get_data(db->memory_handle, index);
        length = reader_get_length(db->memory_handle, index);
    }
    if (db->decompress) {
        std::string pdbText;
        std::string name;
        int err = decompress(data, length, false, pdbText, name);
        if (err != 0) {
            std::string err_msg = "Error decompressing: " + name;
            PyErr_SetString(FoldcompError, err_msg.c_str());
            return NULL;
        }
        PyObject* pdb = PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, pdbText.c_str(), pdbText.size());
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
    (destructor)FoldcompDatabase_dealloc, /* tp_dealloc */
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
    releaseFoldcompDatabase(db);
    Py_RETURN_NONE;
}

static void FoldcompDatabase_dealloc(PyObject* self) {
    FoldcompDatabaseObject* db = (FoldcompDatabaseObject*)self;
    releaseFoldcompDatabase(db);
    Py_TYPE(self)->tp_free(self);
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

// Decompress
int decompress(const char* input, size_t input_size, bool use_alt_order, std::string& output, std::string& name) {
    return decompress(input, input_size, use_alt_order, "pdb", output, name);
}

int decompress(
    const char* input, size_t input_size, bool use_alt_order, const std::string& format,
    std::string& output, std::string& name
) {
    CoutStateGuard coutStateGuard;
    std::cout.setstate(std::ios_base::failbit);
#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
    if (format == "mmcif" || format == "cif") {
        if (!decodeStructureToMMCIF(input, input_size, use_alt_order, name, output)) {
            return 1;
        }
        return 0;
    }
#endif
    if (format != "pdb") {
        return 2;
    }
    if (!decodeStructureToPDB(input, input_size, use_alt_order, name, output)) {
        return 1;
    }
    return 0;
}
// Python binding for decompress
static PyObject *foldcomp_decompress(PyObject* /* self */, PyObject *args, PyObject* kwargs) {
    const char *strArg;
    Py_ssize_t strSize;
    const char* format = "pdb";
    static const char* kwlist[] = {"input", "format", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "y#|$s", const_cast<char**>(kwlist), &strArg, &strSize, &format)) {
        return NULL;
    }

    std::string output;
    std::string name;
    int err = decompress(strArg, strSize, false, format, output, name);
    if (err == 2) {
        PyErr_SetString(PyExc_ValueError, "format must be 'pdb' or 'mmcif'");
        return NULL;
    }
    if (err != 0) {
        PyErr_SetString(FoldcompError, "Error decompressing.");
        return NULL;
    }

    PyObject* pdb = PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, output.c_str(), output.size());
    if (pdb == NULL) {
        return NULL;
    }
    PyObject* result = Py_BuildValue("(s,O)", name.c_str(), pdb);
    Py_DECREF(pdb);
    return result;
}

// Python binding for compress
static PyObject *foldcomp_compress(PyObject* /* self */, PyObject *args, PyObject* kwargs) {
    const char* name;
    const char* pdb_input;
    Py_ssize_t pdb_input_size;
    const char* format = "pdb";
    PyObject* anchor_residue_threshold = NULL;
    PyObject* max_backbone_rmsd = NULL;
    static const char *kwlist[] = {"name", "pdb_content", "format", "anchor_residue_threshold", "max_backbone_rmsd", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss#|$sOO", const_cast<char**>(kwlist),
                                     &name, &pdb_input, &pdb_input_size, &format,
                                     &anchor_residue_threshold, &max_backbone_rmsd)) {
        return NULL;
    }

    if (anchor_residue_threshold != NULL && !PyLong_Check(anchor_residue_threshold)) {
        PyErr_SetString(PyExc_TypeError, "anchor_residue_threshold must be an integer");
        return NULL;
    }
    if (max_backbone_rmsd != NULL &&
        max_backbone_rmsd != Py_None &&
        !PyFloat_Check(max_backbone_rmsd) &&
        !PyLong_Check(max_backbone_rmsd)) {
        PyErr_SetString(PyExc_TypeError, "max_backbone_rmsd must be a float");
        return NULL;
    }

    int threshold = DEFAULT_ANCHOR_THRESHOLD;
    if (anchor_residue_threshold != NULL) {
        threshold = PyLong_AsLong(anchor_residue_threshold);
    }
    float maxBackboneRmsdValue = std::numeric_limits<float>::infinity();
    if (max_backbone_rmsd != NULL && max_backbone_rmsd != Py_None) {
        maxBackboneRmsdValue = static_cast<float>(PyFloat_AsDouble(max_backbone_rmsd));
        if (PyErr_Occurred()) {
            return NULL;
        }
    }

    std::string output;
    int flag = encodeStructureToFoldcompContainer(
        name, pdb_input, static_cast<size_t>(pdb_input_size), format, threshold, maxBackboneRmsdValue, output
    );
    if (flag == PARSE_PDB_NO_ATOM) {
        PyErr_SetString(FoldcompError, "No protein atoms found");
        return NULL;
    } else if (flag == PARSE_PDB_INVALID_FORMAT) {
        PyErr_SetString(PyExc_ValueError, "Invalid structure input or format");
        return NULL;
    } else if (flag != 0) {
        PyErr_SetString(FoldcompError, "Error compressing");
        return NULL;
    }

    return PyBytes_FromStringAndSize(output.c_str(), output.length());
}


PyTypeObject* pathType = NULL;

static PyObject *foldcomp_open(PyObject* /* self */, PyObject* args, PyObject* kwargs) {
    PyObject* path;
    PyObject* user_ids = NULL;
    PyObject* decompress = NULL;
    PyObject* err_on_missing = NULL; // Raise an error if the file is missing. Default: False

    static const char *kwlist[] = {"path", "ids", "decompress", "err_on_missing", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&|$OOO", const_cast<char**>(kwlist), PyUnicode_FSConverter, &path, &user_ids, &decompress, &err_on_missing)) {
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

    std::string dbname(pathCStr);
    std::string index = dbname + ".index";
    bool err_on_missing_flag = false;
    Py_DECREF(path);

    FoldcompDatabaseObject *obj = PyObject_New(FoldcompDatabaseObject, &FoldcompDatabaseType);
    if (obj == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for FoldcompDatabaseObject");
        return NULL;
    }
    obj->user_indices = NULL;
    obj->memory_handle = NULL;
    obj->decompress = true;

    int mode = DB_READER_USE_DATA;
    if (user_ids != NULL && PySequence_Length(user_ids) > 0) {
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

    obj->memory_handle = make_reader(dbname.c_str(), index.c_str(), mode);
    if (obj->memory_handle == NULL) {
        Py_DECREF((PyObject*)obj);
        PyErr_SetString(FoldcompError, "Could not open Foldcomp database");
        return NULL;
    }

    if (user_ids != NULL && PySequence_Length(user_ids) > 0) {
        size_t id_count = (size_t)PySequence_Length(user_ids);
        // Reserve memory for the user indices
        obj->user_indices = new std::vector<int64_t>();
        obj->user_indices->reserve(id_count);
        // user_indices.reserve(id_count);
        for (Py_ssize_t i = 0; i < (Py_ssize_t)id_count; i++) {
            // Iterate over all entries in the database and store ids in a vector of int64_t
            PyObject* item = PySequence_GetItem(user_ids, i);
            if (item == NULL) {
                Py_DECREF((PyObject*)obj);
                return NULL;
            }
            if (!PyUnicode_Check(item)) {
                Py_DECREF(item);
                Py_DECREF((PyObject*)obj);
                PyErr_SetString(PyExc_TypeError, "ids must contain only strings");
                return NULL;
            }
            const char* data = PyUnicode_AsUTF8(item);
            if (data == NULL) {
                Py_DECREF(item);
                Py_DECREF((PyObject*)obj);
                return NULL;
            }
            Py_DECREF(item);
            uint32_t key = reader_lookup_entry(obj->memory_handle, data);
            int64_t id = reader_get_id(obj->memory_handle, key);
            if (id == -1 || key == UINT32_MAX) {
                // Not found --> no error just
                std::string err_msg = "Skipping entry ";
                err_msg += data;
                err_msg += " which is not in the database.";
                if (err_on_missing_flag) {
                    Py_DECREF((PyObject*)obj);
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
    std::string title;
    std::vector<AtomCoordinate> atomCoordinates;
    if (!decodeStructureToAtoms(input, input_size, false, title, atomCoordinates)) {
        PyErr_SetString(PyExc_ValueError, "Could not decompress FCZ file");
        return NULL;
    }

    Foldcomp compRes;
    compRes.compress(atomCoordinates);
    std::vector<float3d> coordsVector = extractCoordinates(atomCoordinates);
    return getPyDictFromFoldcomp(&compRes, coordsVector);
}

// 02. Extract information starting from PDB
PyObject* getDataFromStructureText(const std::string& pdb_input, const char* format) {
    std::vector<AtomCoordinate> atomCoordinates;
    int status = 0;
    if (!parseStructureAtoms(pdb_input.data(), pdb_input.size(), false, atomCoordinates, status, nullptr, format)) {
        if (status == PARSE_PDB_NO_ATOM) {
            PyErr_SetString(PyExc_ValueError, "No protein atoms found in structure input");
            return NULL;
        } else if (status == PARSE_PDB_INVALID_FORMAT) {
            PyErr_SetString(PyExc_ValueError, "Invalid structure input or format");
            return NULL;
        }
        PyErr_SetString(PyExc_ValueError, "Could not parse structure input");
        return NULL;
    }

    Foldcomp compRes;
    compRes.compress(atomCoordinates);

    std::vector<float3d> coordsVector = extractCoordinates(atomCoordinates);

    PyObject* dict = getPyDictFromFoldcomp(&compRes, coordsVector);
    if (dict == NULL) {
        return NULL;
    }
    return dict;
}

static PyObject* foldcomp_get_data(PyObject* /* self */, PyObject* args, PyObject* kwargs) {
    const char* input;
    Py_ssize_t input_size;
    const char* format = "pdb";
    static const char* kwlist[] = { "input", "format", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "y#|$s", (char**)kwlist, &input, &input_size, &format)) {
        return NULL;
    }
    // Check input
    if (input_size == 0) {
        PyErr_SetString(PyExc_ValueError, "Input is empty");
        return NULL;
    }
    // FCMP is standalone FCZ, FCZC is the container format.
    if ((input_size >= MAGICNUMBER_LENGTH && memcmp(input, MAGICNUMBER, MAGICNUMBER_LENGTH) == 0) ||
        hasContainerMagic(input, input_size)) {
        return getDataFromFCZ(input, input_size);
    } else if (input_size >= 4) {
        std::string pdb_input(input, input_size);
        return getDataFromStructureText(pdb_input, format);
    } else {
        PyErr_SetString(PyExc_ValueError, "Input is not a FCZ file or PDB file");
        return NULL;
    }
}

static PyObject* foldcomp_cuda_available(PyObject* /* self */, PyObject* /* args */) {
#ifdef FOLDCOMP_WITH_CUDA
    if (cudaRuntimeAvailable()) Py_RETURN_TRUE;
#endif
    Py_RETURN_FALSE;
}

namespace {

struct PythonBatchResult {
    std::string title;
    std::string text;
};

PyObject* pythonBatchResults(const std::vector<PythonBatchResult>& results) {
    PyObject* list = PyList_New(static_cast<Py_ssize_t>(results.size()));
    if (!list) return nullptr;
    for (size_t i = 0; i < results.size(); ++i) {
        PyObject* text = PyUnicode_FromStringAndSize(
            results[i].text.data(), static_cast<Py_ssize_t>(results[i].text.size()));
        if (!text) {
            Py_DECREF(list);
            return nullptr;
        }
        PyObject* tuple = Py_BuildValue("(s,N)", results[i].title.c_str(), text);
        if (!tuple) {
            Py_DECREF(list);
            return nullptr;
        }
        PyList_SET_ITEM(list, static_cast<Py_ssize_t>(i), tuple);
    }
    return list;
}

} // namespace

#ifdef FOLDCOMP_WITH_CUDA
static PyObject* foldcomp_gpu_pipeline_release(PyObject*, PyObject*) {
    std::string error;
    Py_BEGIN_ALLOW_THREADS
    try {
        releaseGPUPipelineCache();
    } catch (const std::exception& exception) {
        error = exception.what();
    } catch (...) {
        error = "GPU pipeline cache release failed";
    }
    Py_END_ALLOW_THREADS
    if (!error.empty()) {
        PyErr_SetString(FoldcompError, error.c_str());
        return nullptr;
    }
    Py_RETURN_NONE;
}
#endif // FOLDCOMP_WITH_CUDA

static PyObject* foldcomp_decompress_batch(
    PyObject* /* self */, PyObject* args, PyObject* kwargs
) {
    PyObject* inputs;
    PyObject* use_gpu_object = Py_None;
    const char* format = "pdb";
    int max_batch_size = 100;
    int max_residues = 3000;
    int gpu_format = 0;
    int pool_size = 2;
    static const char* kwlist[] = {
        "inputs", "use_gpu", "format", "max_batch_size", "max_residues",
        "gpu_format", "pool_size", nullptr
    };
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O|$Osiipi", const_cast<char**>(kwlist),
            &inputs, &use_gpu_object, &format, &max_batch_size, &max_residues,
            &gpu_format, &pool_size)) {
        return nullptr;
    }
    if (!PyList_Check(inputs)) {
        PyErr_SetString(PyExc_TypeError, "inputs must be a list of bytes objects");
        return nullptr;
    }
    const std::string output_format(format);
    if (output_format != "pdb" && output_format != "mmcif" && output_format != "cif") {
        PyErr_SetString(PyExc_ValueError, "format must be 'pdb' or 'mmcif'");
        return nullptr;
    }
#ifdef FOLDCOMP_WITH_CUDA
    // Not validated here: max_batch_size/max_residues/pool_size only constrain the
    // GPU pipeline, so validating them before use_gpu is resolved below would reject
    // a use_gpu=False call over limits the CPU path (a plain serial loop) never uses.
    // Validated further down, once use_gpu is known to be true.
    GPUDecompressionConfig config;
    config.reader_threads = 1;
    config.max_batch_size = max_batch_size;
    config.max_residues = max_residues;
    config.pool_size = pool_size;
#else
    if (max_batch_size < 1 || max_residues < 1) {
        PyErr_SetString(PyExc_ValueError, "max_batch_size and max_residues must be positive");
        return nullptr;
    }
#endif

    const Py_ssize_t count = PyList_GET_SIZE(inputs);
    std::vector<std::string> buffers;
    buffers.reserve(static_cast<size_t>(count));
    for (Py_ssize_t i = 0; i < count; ++i) {
        PyObject* item = PyList_GET_ITEM(inputs, i);
        if (!PyBytes_Check(item)) {
            PyErr_Format(PyExc_TypeError, "inputs[%zd] must be bytes", i);
            return nullptr;
        }
        buffers.emplace_back(
            PyBytes_AS_STRING(item), static_cast<size_t>(PyBytes_GET_SIZE(item)));
    }

    bool runtime_available = false;
#ifdef FOLDCOMP_WITH_CUDA
    runtime_available = cudaRuntimeAvailable();
#endif
    bool use_gpu = runtime_available;
    if (use_gpu_object != Py_None) {
        if (!PyBool_Check(use_gpu_object)) {
            PyErr_SetString(PyExc_TypeError, "use_gpu must be bool or None");
            return nullptr;
        }
        int truth = PyObject_IsTrue(use_gpu_object);
        if (truth < 0) return nullptr;
        use_gpu = truth != 0;
    }
    if (use_gpu && !runtime_available) {
        PyErr_SetString(FoldcompError, "CUDA GPU support is not available at runtime");
        return nullptr;
    }
#ifdef FOLDCOMP_WITH_CUDA
    if (use_gpu) {
        std::string config_error;
        if (!config.validate(config_error)) {
            PyErr_SetString(PyExc_ValueError, config_error.c_str());
            return nullptr;
        }
    }
#endif

    std::vector<PythonBatchResult> results(static_cast<size_t>(count));
    std::string error;
    Py_ssize_t failed_index = -1;

    if (!use_gpu) {
        Py_BEGIN_ALLOW_THREADS
        try {
            for (Py_ssize_t i = 0; i < count; ++i) {
                if (decompress(buffers[i].data(), buffers[i].size(), false, output_format,
                               results[i].text, results[i].title) != 0) {
                    failed_index = i;
                    break;
                }
            }
        } catch (const std::exception& exception) {
            error = exception.what();
        } catch (...) {
            error = "CPU batch decompression failed";
        }
        Py_END_ALLOW_THREADS
    }
#ifdef FOLDCOMP_WITH_CUDA
    else {
        std::vector<MemoryViewProcessor::Entry> entries;
        entries.reserve(static_cast<size_t>(count));
        for (Py_ssize_t i = 0; i < count; ++i) {
            entries.emplace_back(
                std::to_string(i),
                tcb::span<const char>(buffers[i].data(), buffers[i].size()));
        }
        MemoryViewProcessor processor(entries);
        // config (reader_threads/max_batch_size/max_residues) was already built and
        // validated above, once use_gpu was known to be true (this branch).
        // GPU PDB writer: format ATOM lines on the device in-pipeline (PDB only).
        const bool use_gpu_pdb = gpu_format && output_format == "pdb";
        config.emit_device_pdb = use_gpu_pdb;
        std::vector<bool> completed(static_cast<size_t>(count), false);
        // For the GPU PDB path the per-atom text assembly is the wall-clock
        // bottleneck (the GPU is otherwise ~80% idle), so defer it out of the
        // single drain thread and stitch all inputs in parallel after the run.
        std::vector<GPUDecompressionResult> pdb_staging(
            use_gpu_pdb ? static_cast<size_t>(count) : 0);
        auto collect = [&](GPUDecompressionResult&& result) -> bool {
            size_t index = 0;
            try {
                index = static_cast<size_t>(std::stoul(result.name));
            } catch (...) {
                return false;
            }
            if (index >= results.size()) return false;
            results[index].title = result.title;
#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
            if (output_format == "mmcif" || output_format == "cif") {
                writeSegmentsToMMCIF(result.segments, result.title, results[index].text);
            } else
#endif
            if (use_gpu_pdb) {
                // Defer: just capture the GPU-formatted lines + metadata (cheap
                // move); the stitch runs in parallel after the run.
                pdb_staging[index] = std::move(result);
            } else {
                writeSegmentsToPDB(result.segments, result.title, results[index].text);
            }
            completed[index] = true;
            return true;
        };
        bool success = false;
        Py_BEGIN_ALLOW_THREADS
        try {
            GPUPipelineLease lease; // scoped to the try block: returns to the pool as
                                     // soon as the run + PDB stitch below are done.
            {
                nvtx3::scoped_range r{"GPU pipeline (run)"};
                success = runGPUDecompressionPipeline(processor, config, collect, error, &lease);
            }
            if (success && use_gpu_pdb) {
                // Parallel PDB assembly: each input is an independent string read
                // straight from the GPU line-capture buffer of this run's pipeline
                // instance (no serial drain copy). Byte-identical to the serial
                // version; GIL released here. The parallel stitch is a large,
                // size-independent win over serial (measured ~3ms flat vs serial's
                // linear growth to ~28ms at N=256), so it always runs across the
                // team. The idle workers' post-region GOMP keep-alive spin overlaps
                // the subsequent single-threaded pythonBatchResults build and does
                // not add wall-clock.
                nvtx3::scoped_range r{"GPU PDB stitch (host)"};
                // `lease` keeps this run's pipeline instance checked out of the
                // pool for as long as it's alive, so no other thread can reuse/evict
                // it (and free its capture buffer) while this stitch reads from it.
                const char* base = lease.pdbCaptureBase();
                const Py_ssize_t n = count;
                // An exception escaping an OpenMP structured block is undefined
                // behavior (std::terminate), so a bad_alloc or other failure inside
                // the loop body must be caught per-iteration and re-thrown only
                // after the region ends, from single-threaded code where the
                // enclosing try/catch below can turn it into `error`.
                std::atomic<bool> stitch_failed{false};
                std::mutex stitch_error_mutex;
                std::string stitch_error;
#pragma omp parallel for schedule(dynamic)
                for (Py_ssize_t i = 0; i < n; ++i) {
                    if (stitch_failed.load(std::memory_order_relaxed)) {
                        continue;
                    }
                    try {
                        const GPUDecompressionResult& r = pdb_staging[static_cast<size_t>(i)];
                        std::string& out = results[static_cast<size_t>(i)].text;
                        if (r.pdb_frags.size() == 1) {
                            // Single-fragment input: lines are contiguous in the buffer.
                            writeLineSegmentsToPDB(
                                r.line_segments, r.title,
                                base + r.pdb_frags[0].offset * PDB_ATOM_LINE_LEN, out);
                        } else {
                            // Multi-fragment container: gather the (non-contiguous)
                            // fragments into a temporary, then stitch.
                            size_t tot = 0;
                            for (const GPUCoordFragment& f : r.pdb_frags)
                                tot += static_cast<size_t>(f.n_atoms);
                            std::string tmp(tot * PDB_ATOM_LINE_LEN, '\0');
                            size_t o = 0;
                            for (const GPUCoordFragment& f : r.pdb_frags) {
                                std::memcpy(&tmp[o * PDB_ATOM_LINE_LEN],
                                    base + f.offset * PDB_ATOM_LINE_LEN,
                                    static_cast<size_t>(f.n_atoms) * PDB_ATOM_LINE_LEN);
                                o += static_cast<size_t>(f.n_atoms);
                            }
                            writeLineSegmentsToPDB(r.line_segments, r.title, tmp.data(), out);
                        }
                    } catch (const std::exception& exception) {
                        if (!stitch_failed.exchange(true)) {
                            std::lock_guard<std::mutex> lock(stitch_error_mutex);
                            stitch_error = exception.what();
                        }
                    } catch (...) {
                        if (!stitch_failed.exchange(true)) {
                            std::lock_guard<std::mutex> lock(stitch_error_mutex);
                            stitch_error = "GPU PDB stitch failed";
                        }
                    }
                }
                if (stitch_failed.load()) {
                    throw std::runtime_error(stitch_error);
                }
            }
        } catch (const std::exception& exception) {
            error = exception.what();
        } catch (...) {
            error = "GPU batch decompression failed";
        }
        Py_END_ALLOW_THREADS
        if (success) {
            for (size_t i = 0; i < completed.size(); ++i) {
                if (!completed[i]) {
                    success = false;
                    error = "GPU pipeline did not emit input " + std::to_string(i);
                    break;
                }
            }
        }
        if (!success && error.empty()) error = "GPU batch decompression failed";
    }
#endif

    if (failed_index >= 0) {
        PyErr_Format(FoldcompError, "Error decompressing input %zd", failed_index);
        return nullptr;
    }
    if (!error.empty()) {
        PyErr_SetString(FoldcompError, error.c_str());
        return nullptr;
    }
    FOLDCOMP_CPU_NVTX(nvtx_build, "pythonBatchResults (PyUnicode)");
    return pythonBatchResults(results);
}

// Method definitions
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wcast-function-type"
static PyMethodDef foldcomp_methods[] = {
    // {"compress", foldcomp_compress, METH_VARARGS, "Compress a PDB file."},
    {"decompress", (PyCFunction)foldcomp_decompress, METH_VARARGS | METH_KEYWORDS, "Decompress Foldcomp content to PDB or mmCIF."},
    {"compress", (PyCFunction)foldcomp_compress, METH_VARARGS | METH_KEYWORDS, "Compress PDB content to FCZ."},
    {"open", (PyCFunction)foldcomp_open, METH_VARARGS | METH_KEYWORDS, "Open a Foldcomp database."},
    {"get_data", (PyCFunction)foldcomp_get_data, METH_VARARGS | METH_KEYWORDS, "Get data from FCZ or PDB content."},
    {"cuda_available", foldcomp_cuda_available, METH_NOARGS, "Return whether a usable CUDA device is available."},
    {"decompress_batch", (PyCFunction)foldcomp_decompress_batch,
     METH_VARARGS | METH_KEYWORDS, "Batch-decompress FCMP/FCZC byte strings."},
#ifdef FOLDCOMP_WITH_CUDA
    {"gpu_pipeline_release", foldcomp_gpu_pipeline_release, METH_NOARGS,
     "Release the cached GPU decompression pipeline (persistent slots + pinned "
     "buffers reused across decompress_batch(use_gpu=True) calls). The next "
     "call rebuilds it on demand."},
#endif
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

    // Keep the module-specific exception while making it compatible with the
    // standard exception type promised by the batch API.
    FoldcompError = PyErr_NewException("foldcomp.error", PyExc_RuntimeError, NULL);
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
