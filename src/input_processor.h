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
#include <vector>
#include <string>
#include <iterator>
#include <atomic>
#include <mutex>
#include <cstdlib>

#ifdef HAVE_GCS
#include "google/cloud/storage/client.h"
namespace gcs = ::google::cloud::storage;
#endif

// #ifdef HAVE_AWS_S3
// #include <aws/core/Aws.h>
// #include <aws/s3/S3Client.h>
// #include <aws/s3/model/ListObjectsRequest.h>
// #include <aws/s3/model/GetObjectRequest.h>
// #endif

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
    DirectoryProcessor(std::vector<std::string> files) : files(std::move(files)) {};

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
        int mode = DB_READER_USE_DATA | DB_READER_USE_LOOKUP_REVERSE;
        handle = make_reader(input.c_str(), index.c_str(), mode);
        id_as_name = 1;
    };
    
    DatabaseProcessor(const std::string& input, std::string& user_id_file, int id_mode, int use_cache) {
        std::string index = input + ".index";
        id_as_name = id_mode;
        int mode = DB_READER_USE_DATA;
        if (use_cache) {
            mode |= DB_READER_CACHE;
        } else if (id_as_name) {
            mode |= DB_READER_USE_LOOKUP;
        }
        
        handle = make_reader(input.c_str(), index.c_str(), mode);
        _read_id_list(user_id_file);
    };

    DatabaseProcessor(const std::string& input, const std::vector<std::string>& ids, int id_mode) {
        std::string index = input + ".index";
        int mode = DB_READER_USE_DATA | DB_READER_USE_LOOKUP;
        handle = make_reader(input.c_str(), index.c_str(), mode);
        user_ids = ids;
        id_as_name = id_mode;
    };

    ~DatabaseProcessor() {
        free_reader(handle);
    };

    void run(process_entry_func func, int num_threads) override {
        size_t db_size = reader_get_size(handle);
#pragma omp parallel shared(handle) num_threads(num_threads)
        if (user_ids.size() == 0) { // process all entries in db
            {
#pragma omp for
                for (size_t i = 0; i < db_size; i++) {
                    uint32_t key = reader_get_key(handle, i);
                    const char* name = reader_lookup_name_alloc(handle, key);
                    // If name == "", throw warning
                    if (name && !name[0]) {
                        std::cerr << "[Warning] empty name for key " << key << std::endl;
                        free((void*)name);
                        continue;
                    }
                    if (!func(name, reader_get_data(handle, i), reader_get_length(handle, i))) {
                        std::cerr << "[Error] processing db entry " << name << " failed." << std::endl;
                        free((void*)name);
                        continue;
                    }
                    free((void*)name);
                }
            }
        } else { // process only entries in user_ids
#pragma omp for
            for (size_t i = 0; i < user_ids.size(); i++) {
                uint32_t key;

                if (!id_as_name) {
                    key = std::stoi(user_ids[i]);
                } else {
                    key = reader_lookup_entry(handle, user_ids[i].c_str());
                }
                int64_t id = reader_get_id(handle, key);
                if (id == -1 || key == UINT32_MAX) {
                    // NOT found
                    std::cerr << "[Warning] " << user_ids[i] << " not found in database." << std::endl;
                    continue;
                }
                if (!func(user_ids[i].c_str(), reader_get_data(handle, id), reader_get_length(handle, id))) {
                    std::cerr << "[Error] processing db entry " << user_ids[i] << " failed." << std::endl;
                    continue;
                }
            }
        }
    }

private:
    void* handle;
    std::vector<std::string> user_ids;
    int id_as_name;

    void _read_id_list(std::string& file) {
        // Check if file exists
        if (!std::ifstream(file)) {
            std::cerr << "[Error] user id '" << file << "' does not exist." << std::endl;
            return;
        }
        std::ifstream infile(file);
        std::string line;
        while (std::getline(infile, line)) {
            user_ids.push_back(line);
        }
        infile.close();
    }
};

#ifdef HAVE_GCS
class GcsProcessor : public Processor {
public:
    GcsProcessor(std::vector<std::string> object_uris) : object_uris(std::move(object_uris)) {};
    GcsProcessor(const std::string& object_uri) : GcsProcessor(std::vector<std::string>{ object_uri }) {};

    void run(process_entry_func func, int num_threads) override {
        const int worker_threads = num_threads > 0 ? num_threads : 1;
        const char* adc_path = std::getenv("GOOGLE_APPLICATION_CREDENTIALS");

        if (adc_path == NULL || adc_path[0] == '\0') {
            log_message("GOOGLE_APPLICATION_CREDENTIALS is not set. "
                        "If startup appears to hang, set ADC explicitly to your service-account JSON.");
        }

        // Log volume scales with the input instead of being pinned to a fixed cadence:
        // `block` is the largest power of ten <= total/100 and ends a line with the
        // running count, and 25 '=' fill that block, so one '=' is block/25 objects.
        // A '=' never stands for fewer than 10 objects; when that floor binds (inputs
        // under ~25k) the block follows it rather than the other way round.
        const size_t total = object_uris.size();
        size_t block = 1;
        while (block * 10 <= total / 100) block *= 10;
        size_t tick = block / 25;
        if (tick < 10) {
            tick = 10;
            block = 25 * tick;
        }

        log_message("Starting GCS processing for " + std::to_string(total) +
                    " object(s) with " + std::to_string(worker_threads) + " thread(s)");

        auto options = google::cloud::Options{}
            .set<gcs::ConnectionPoolSizeOption>(worker_threads)
            .set<google::cloud::storage_experimental::HttpVersionOption>("2.0");
        gcs::Client client(options);
        log_message("GCS client initialized");

        std::atomic<size_t> processed_count(0);
        std::atomic<size_t> success_count(0);
        std::atomic<size_t> failed_count(0);

#pragma omp parallel for schedule(dynamic) num_threads(worker_threads)
        for (size_t i = 0; i < object_uris.size(); ++i) {
            bool success = false;
            std::string bucket_name;
            std::string object_name;
            if (!parse_uri(object_uris[i], bucket_name, object_name)) {
                std::cerr << "[Error] Invalid GCS URI: " << object_uris[i] << std::endl;
            } else {
                auto reader = client.ReadObject(bucket_name, object_name);
                if (!reader.status().ok()) {
                    std::cerr << "[Error] Could not read GCS object " << object_uris[i]
                              << ": " << reader.status() << std::endl;
                } else {
                    std::string contents((std::istreambuf_iterator<char>(reader)), std::istreambuf_iterator<char>());
                    success = func(object_name.c_str(), contents.c_str(), contents.size());
                    if (!success) {
                        std::cerr << "[Error] processing GCS object " << object_uris[i] << " failed." << std::endl;
                    }
                }
            }

            if (success) {
                success_count.fetch_add(1, std::memory_order_relaxed);
            } else {
                failed_count.fetch_add(1, std::memory_order_relaxed);
            }

            // Same shape as mmseqs' Debug::Progress on a non-tty: a run of '=' capped by
            // a tab and the running count. The boundary writes its own '=' before the
            // label, so a full line carries exactly 25.
            size_t done = processed_count.fetch_add(1, std::memory_order_relaxed) + 1;
            if (done % block == 0) {
                progress_write("=\t" + human_count(done) + " structures processed\n");
            } else if (done % tick == 0) {
                progress_write("=");
            }
        }

        // Close whatever partial '=' run the last block left open, so the summary does
        // not get appended to it.
        if (total >= tick && total % block != 0) {
            progress_write("\n");
        }
        log_message(
            "Finished GCS processing: ok=" +
            std::to_string(success_count.load(std::memory_order_relaxed)) +
            ", failed=" +
            std::to_string(failed_count.load(std::memory_order_relaxed))
        );
    }

private:
    static std::mutex& log_mutex() {
        static std::mutex m;
        return m;
    }

    static void log_message(const std::string& message) {
        std::lock_guard<std::mutex> lock(log_mutex());
        std::cerr << "[GCS] " << message << std::endl;
    }

    // Progress is written raw — no "[GCS] " prefix and no newline — because the bar is
    // built up one '=' at a time across a line. Shares log_message's mutex so the two
    // cannot interleave mid-line.
    static void progress_write(const std::string& s) {
        std::lock_guard<std::mutex> lock(log_mutex());
        std::cerr << s;
        std::cerr.flush();
    }

    // Counts land on a power of ten, so an exact suffix is always available.
    static std::string human_count(size_t n) {
        if (n >= 1000000 && n % 1000000 == 0) return std::to_string(n / 1000000) + " Mio.";
        if (n >= 1000 && n % 1000 == 0) return std::to_string(n / 1000) + " K";
        return std::to_string(n);
    }

    static bool parse_uri(const std::string& uri, std::string& bucket_name, std::string& object_name) {
        std::string path;
        if (stringStartsWith("gcs://", uri)) {
            path = uri.substr(std::string("gcs://").size());
        } else if (stringStartsWith("gs://", uri)) {
            path = uri.substr(std::string("gs://").size());
        } else {
            return false;
        }
        const size_t slash = path.find('/');
        if (slash == std::string::npos || slash == 0 || slash + 1 == path.size()) {
            return false;
        }
        bucket_name = path.substr(0, slash);
        object_name = path.substr(slash + 1);
        return true;
    }

    std::vector<std::string> object_uris;
};
#endif

// #ifdef HAVE_AWS_S3
// class S3Processor : public Processor {
// public:
//     S3Processor(const std::string& input) {
//         Aws::Client::ClientConfiguration clientConfig;
//         clientConfig.region = Aws::Region::US_WEST_2;  // Update the region if needed
//         client = std::make_shared<Aws::S3::S3Client>(clientConfig);
//         bucket_name = input;
//     };

//     void run(process_entry_func func, int num_threads) override {
// #pragma omp parallel num_threads(num_threads)
//         {
// #pragma omp single
//             // Get object list from S3 bucket
//             Aws::S3::Model::ListObjectsRequest objects_request;
//             objects_request.WithBucket(bucket_name);

//             auto list_objects_outcome = client->ListObjects(objects_request);

//             if (list_objects_outcome.IsSuccess()) {
//                 auto object_list = list_objects_outcome.GetResult().GetContents();
//                 for (auto const& s3_object : object_list) {
//                     std::string obj_name = s3_object.GetKey();
// #pragma omp task firstprivate(obj_name)
//                     {
//                         bool skipFilter = true;
//                         bool allowedSuffix = stringEndsWith(".cif", obj_name) || stringEndsWith(".pdb", obj_name);
//                         if (skipFilter && allowedSuffix) {
//                             Aws::S3::Model::GetObjectRequest object_request;
//                             object_request.WithBucket(bucket_name).WithKey(obj_name);

//                             auto get_object_outcome = client->GetObject(object_request);

//                             if (get_object_outcome.IsSuccess()) {
//                                 auto& retrieved_file = get_object_outcome.GetResultWithOwnership().GetBody();
//                                 std::string contents{ std::istreambuf_iterator<char>{retrieved_file}, {} };
//                                 func(obj_name.c_str(), contents.c_str(), contents.length());
//                             }
//                             else {
//                                 std::cerr << "Could not read object " << obj_name << std::endl;
//                             }
//                         }
//                     }
//                 }
//             }
//             else {
//                 std::cout << "ListObjects error: "
//                     << list_objects_outcome.GetError().GetExceptionName() << " - "
//                     << list_objects_outcome.GetError().GetMessage() << std::endl;
//             }
//         }
//     }

// private:
//     std::shared_ptr<Aws::S3::S3Client> client;
//     std::string bucket_name;
// };
// #endif
