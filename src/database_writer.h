/**
 * File: database_writer.h
 * Created: 2022-12-09 14:53:33
 * Author: Milot Mirdita (milot@mirdita.de)
 */

#pragma once
#include <cstdint>
#include <cstddef>

void* make_writer(const char *data_name, const char *index_name);
void free_writer(void *reader);

bool writer_append(void *reader, const char* data, size_t length, uint32_t key, const char* name);
