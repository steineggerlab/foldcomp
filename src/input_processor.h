/**
 * File: input_processor.h
 * Created: 2023-02-10 17:04:08
 * Author: Milot Mirdita (milot@mirdita.de)
 */

#pragma once

#include "microtar.h"
#include "utility.h"
#include "database_reader.h"

#include <utility>
#include <functional>

#ifdef HAVE_GCS
#include "google/cloud/storage/client.h"
#endif

// OpenMP for parallelization
#ifdef OPENMP
#include <omp.h>
#endif

#include <zlib.h>
static int file_gzread(mtar_t *tar, void *data, size_t size) {
    size_t res = gzread((gzFile)tar->stream, data, size);
    return (res == size) ? MTAR_ESUCCESS : MTAR_EREADFAIL;
}

static int file_gzseek(mtar_t *tar, long offset, int whence) {
    int res = gzseek((gzFile)tar->stream, offset, whence);
    return (res != -1) ? MTAR_ESUCCESS : MTAR_ESEEKFAIL;
}

static int file_gzclose(mtar_t *tar) {
    gzclose((gzFile)tar->stream);
    return MTAR_ESUCCESS;
}

int mtar_gzopen(mtar_t *tar, const char *filename) {
    // Init tar struct and functions
    memset(tar, 0, sizeof(*tar));
    tar->read = file_gzread;
    tar->seek = file_gzseek;
    tar->close = file_gzclose;
    // Open file
    tar->stream = gzopen(filename, "rb");
    if (!tar->stream) {
        return MTAR_EOPENFAIL;
    }

#if defined(ZLIB_VERNUM) && ZLIB_VERNUM >= 0x1240
    gzbuffer((gzFile)tar->stream, 1 * 1024 * 1024);
#endif

    return MTAR_ESUCCESS;
}

using process_entry_func = std::function<bool(const char* name, const char* content, size_t size)>;

class Processor {
public:
    virtual ~Processor() {};
    virtual void run(process_entry_func, int) {};
};

class DirectoryProcessor : public Processor {
public:
    DirectoryProcessor(const std::string& input, bool recursive) {
        files = getFilesInDirectory(input, recursive);
    };

    void run(process_entry_func func, int num_threads) override {
#pragma omp parallel shared(files) num_threads(num_threads)
        {
            char* dataBuffer;
            ssize_t bufferSize;
#pragma omp for
            for (size_t i = 0; i < files.size(); i++) {
                std::string name = files[i];
                FILE* file = fopen(name.c_str(), "r");
                dataBuffer = file_map(file, &bufferSize, 0);
                if (!func(name.c_str(), dataBuffer, bufferSize)) {
                    std::cerr << "[Error] processing dir entry " << name << " failed." << std::endl;
                    file_unmap(dataBuffer, bufferSize);
                    continue;
                }
                file_unmap(dataBuffer, bufferSize);
                fclose(file);
            }
        }
    }

private:
    std::vector<std::string> files;
};

class TarProcessor : public Processor {
public:
    TarProcessor(const std::string& input) {
        if (stringEndsWith(".gz", input) || stringEndsWith(".tgz", input)) {
            if (mtar_gzopen(&tar, input.c_str()) != MTAR_ESUCCESS) {
                std::cerr << "[Error] open tar " << input << " failed." << std::endl;
            }
        } else {
            if (mtar_open(&tar, input.c_str(), "r") != MTAR_ESUCCESS) {
                std::cerr << "[Error] open tar " << input << " failed." << std::endl;
            }
        }
    };

    ~TarProcessor() {
        mtar_close(&tar);
    };

    void run(process_entry_func func, int num_threads) override {
#pragma omp parallel shared(tar) num_threads(num_threads)
        {
            bool proceed = true;
            mtar_header_t header;
            size_t bufferSize = 1024 * 1024;
            char* dataBuffer = (char*)malloc(bufferSize);
            std::string name;
            while (proceed) {
                bool writeEntry = true;
#pragma omp critical
                {
                    if (tar.isFinished == 0 && (mtar_read_header(&tar, &header)) != MTAR_ENULLRECORD) {
                        // GNU tar has special blocks for long filenames
                        if (header.type == MTAR_TGNU_LONGNAME || header.type == MTAR_TGNU_LONGLINK) {
                            if (header.size > bufferSize) {
                                bufferSize = header.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                std::cerr << "[Error] cannot read entry " << header.name << std::endl;
                                goto done;
                            }
                            name.assign(dataBuffer, header.size);
                            // skip to next record
                            if (mtar_read_header(&tar, &header) == MTAR_ENULLRECORD) {
                                std::cerr << "[Error] tar truncated after entry " << name << std::endl;
                                goto done;
                            }
                        } else {
                            name = header.name;
                        }
                        if (header.type == MTAR_TREG || header.type == MTAR_TCONT || header.type == MTAR_TOLDREG) {
                            if (header.size > bufferSize) {
                                bufferSize = header.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                std::cerr << "[Error] cannot read entry " << name << std::endl;
                                goto done;
                            }
                            proceed = true;
                            writeEntry = true;
                        } else {
                            if (header.size > 0 && mtar_skip_data(&tar) != MTAR_ESUCCESS) {
                                std::cerr << "[Error] cannot skip entry " << name << std::endl;
                                goto done;
                            }
                            proceed = true;
                            writeEntry = false;
                        }
                    } else {
done:
                        tar.isFinished = 1;
                        proceed = false;
                        writeEntry = false;
                    }
                } // end read in
                if (proceed && writeEntry) {
                    if (!func(name.c_str(), dataBuffer, header.size)) {
                        std::cerr << "[Error] failed processing tar entry " << name << std::endl;
                        continue;
                    }
                }
            }
            free(dataBuffer);
        }
    }

private:
    mtar_t tar;
};

class DatabaseProcessor : public Processor {
public:
    DatabaseProcessor(const std::string& input) {
        std::string index = input + ".index";
        int mode = 0;
        handle = make_reader(input.c_str(), index.c_str(), mode);
    };
    ~DatabaseProcessor() {
        free_reader(handle);
    };

    void run(process_entry_func func, int num_threads) override {
        size_t db_size = reader_get_size(handle);
#pragma omp parallel shared(handle) num_threads(num_threads)
        {
#pragma omp for
            for (size_t i = 0; i < db_size; i++) {
                // TODO
                // std::string name = reader_get_name(handle, i);
                std::string name = std::to_string(i);
                std::cout << "processing " << name << std::endl;
                if (!func(name.c_str(), reader_get_data(handle, i), reader_get_length(handle, i))) {
                    std::cerr << "[Error] processing db entry " << name << " failed." << std::endl;
                    continue;
                }
            }
        }
    }

private:
    void* handle;
};

#ifdef HAVE_GCS
class GcsProcessor : public Processor {
public:
    namespace gcs = ::google::cloud::storage;
    GcsProcessor(const std::string& input) {
        auto options = google::cloud::Options{}
            .set<gcs::ConnectionPoolSizeOption>(num_threads)
            .set<google::cloud::storage_experimental::HttpVersionOption>("2.0");
        client = gcs::Client(options);
        bucket_name = input;
    };

    void run(process_entry_func func, int num_threads) override {
       #pragma omp parallel num_threads(num_threads)
        {
#pragma omp single
            // Get object list from gcs bucket
            for (auto&& object_metadata : client.ListObjects(bucket_name, gcs::Projection::NoAcl(), gcs::MaxResults(100000))) {
                std::string obj_name = object_metadata->name();
                // Set zero padding for ID with 4 digits
#pragma omp task firstprivate(obj_name)
                {
                    // Filter for splitting input into 10 different processes
                    // bool skipFilter = filter != '\0' && obj_name.length() >= 9 && obj_name[8] == filter;
                    bool skipFilter = true;
                    bool allowedSuffix = stringEndsWith(".cif", obj_name) || stringEndsWith(".pdb", obj_name);
                    if (skipFilter && allowedSuffix) {
                        auto reader = client.ReadObject(bucket_name, obj_name);
                        if (!reader.status().ok()) {
                            std::cerr << "Could not read object " << obj_name << std::endl;
                        } else {
                            std::string contents{ std::istreambuf_iterator<char>{reader}, {} };
                            func(obj_name.c_str(), contents.c_str(), contents.length());
                        }
                    }
                }
            }
        }
    }

private:
    google::cloud::storage::Client::Client client;
    std::string bucket_name;
};
#endif
