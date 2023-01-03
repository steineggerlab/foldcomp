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
        if (mtar_open(&tar, input.c_str(), "r") != MTAR_ESUCCESS) {
            std::cerr << "[Error] open tar " << input << " failed." << std::endl;
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
                    if (mtar_read_header(&tar, &header) != MTAR_ENULLRECORD) {
                        //TODO GNU tar has special blocks for long filenames
                        name = header.name;
                        if (header.size > bufferSize) {
                            bufferSize = header.size * 1.5;
                            dataBuffer = (char*)realloc(dataBuffer, bufferSize);
                        }
                        if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                            std::cerr << "[Error] reading tar entry " << name << " failed." << std::endl;
                            writeEntry = false;
                            proceed = false;
                        }
                        else {
                            writeEntry = true;
                            proceed = true;
                        }
                        mtar_next(&tar);
                        writeEntry = (header.type == MTAR_TREG) ? writeEntry : false;
                    }
                    else {
                        proceed = false;
                        writeEntry = false;
                    }
                } // end read in
                if (proceed && writeEntry) {
                    if (!func(name.c_str(), dataBuffer, header.size)) {
                        std::cerr << "[Error] processing tar entry " << name << " failed." << std::endl;
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
