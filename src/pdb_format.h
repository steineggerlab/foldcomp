// SPDX-License-Identifier: MIT
/**
 * File: pdb_format.h
 * Project: foldcomp
 * Description:
 *     Cursor-based PDB column formatters usable from both host and device code.
 *     `pf_format_atom_line` writes exactly PDB_ATOM_LINE_LEN (81) bytes and MUST
 *     stay byte-for-byte identical to the ATOM-record formatting in
 *     writeAtomCoordinatesToPDB (src/atom_coordinate.cpp) — the GPU PDB writer is
 *     verified against that host formatter.
 */
#pragma once
#include "atom_coordinate.h"
#include <cstddef>
#include <cstdint>

#if defined(__CUDACC__)
#define PDB_HD __host__ __device__
#else
#define PDB_HD
#endif

// PDB_ATOM_LINE_LEN is defined in atom_coordinate.h (included above).

PDB_HD inline size_t pf_strlen(const char* s) {
    size_t n = 0;
    while (s[n]) ++n;
    return n;
}

PDB_HD inline char* pf_left(char* p, const char* v, size_t vlen, size_t width) {
    size_t copied = vlen < width ? vlen : width;
    for (size_t i = 0; i < copied; ++i) *p++ = v[i];
    for (size_t i = copied; i < width; ++i) *p++ = ' ';
    return p;
}

PDB_HD inline char* pf_right_str(char* p, const char* v, size_t vlen, size_t width) {
    size_t copied = vlen < width ? vlen : width;
    for (size_t i = copied; i < width; ++i) *p++ = ' ';
    for (size_t i = 0; i < copied; ++i) *p++ = v[vlen - copied + i];
    return p;
}

PDB_HD inline char* pf_int(char* p, int value, size_t width) {
    bool negative = value < 0;
    uint32_t magnitude = negative
        ? static_cast<uint32_t>(-(static_cast<int64_t>(value)))
        : static_cast<uint32_t>(value);
    char digits[16];
    int len = 0;
    do {
        digits[len++] = static_cast<char>('0' + (magnitude % 10));
        magnitude /= 10;
    } while (magnitude != 0);
    size_t totalLen = static_cast<size_t>(len + (negative ? 1 : 0));
    if (totalLen > width) {
        // Value doesn't fit in the fixed-width field: emit `width` overflow
        // markers instead of writing past the caller's fixed-size record.
        for (size_t i = 0; i < width; ++i) *p++ = '*';
        return p;
    }
    for (size_t i = totalLen; i < width; ++i) *p++ = ' ';
    if (negative) *p++ = '-';
    while (len-- > 0) *p++ = digits[len];
    return p;
}

template <size_t Width, uint32_t Scale, size_t Precision>
PDB_HD inline char* pf_fixed(char* p, float value) {
    bool negative = value < 0.0f;
    float adjusted = value + (negative ? -(0.5f / static_cast<float>(Scale))
                                       :  (0.5f / static_cast<float>(Scale)));
    int64_t scaled = static_cast<int64_t>(adjusted * static_cast<float>(Scale));
    if (scaled < 0) scaled = -scaled;
    uint32_t fraction = static_cast<uint32_t>(scaled % Scale);
    uint32_t integer = static_cast<uint32_t>(scaled / Scale);

    char intDigits[16];
    int intLen = 0;
    do {
        intDigits[intLen++] = static_cast<char>('0' + (integer % 10));
        integer /= 10;
    } while (integer != 0);

    size_t totalLen = static_cast<size_t>(intLen) + 1 + Precision + (negative ? 1 : 0);
    if (totalLen > Width) {
        // Value doesn't fit in the fixed-width field: emit `Width` overflow
        // markers instead of writing past the caller's fixed-size record.
        for (size_t i = 0; i < Width; ++i) *p++ = '*';
        return p;
    }
    for (size_t i = totalLen; i < Width; ++i) *p++ = ' ';
    if (negative) *p++ = '-';
    while (intLen-- > 0) *p++ = intDigits[intLen];
    *p++ = '.';

    char fracDigits[Precision];
    for (size_t i = 0; i < Precision; ++i) {
        fracDigits[Precision - 1 - i] = static_cast<char>('0' + (fraction % 10));
        fraction /= 10;
    }
    for (size_t i = 0; i < Precision; ++i) *p++ = fracDigits[i];
    return p;
}

/**
 * Write one PDB ATOM record (exactly PDB_ATOM_LINE_LEN bytes, trailing '\n') to
 * `dst`. Mirrors the ATOM-record block of writeAtomCoordinatesToPDB.
 */
PDB_HD inline void pf_format_atom_line(const AtomCoordinate& a, char* dst) {
    char* p = dst;
    // "ATOM  "
    *p++ = 'A'; *p++ = 'T'; *p++ = 'O'; *p++ = 'M'; *p++ = ' '; *p++ = ' ';
    p = pf_int(p, a.atom_index, 5);
    *p++ = ' ';
    size_t nameLen = pf_strlen(a.atom.data);
    if (nameLen == 4) {
        p = pf_left(p, a.atom.data, nameLen, 4);
    } else {
        *p++ = ' ';
        p = pf_left(p, a.atom.data, nameLen, 3);
    }
    *p++ = (a.altloc == '\0') ? ' ' : a.altloc;
    p = pf_right_str(p, a.residue.data, pf_strlen(a.residue.data), 3);
    *p++ = ' ';
    char chainId = (a.chain.data[0] == '\0') ? ' ' : a.chain.data[0];
    *p++ = chainId;
    p = pf_int(p, a.residue_index, 4);
    *p++ = (a.insertion_code == '\0') ? ' ' : a.insertion_code;
    *p++ = ' '; *p++ = ' '; *p++ = ' ';
    p = pf_fixed<8, 1000, 3>(p, a.coordinate.x);
    p = pf_fixed<8, 1000, 3>(p, a.coordinate.y);
    p = pf_fixed<8, 1000, 3>(p, a.coordinate.z);
    p = pf_fixed<6, 100, 2>(p, a.occupancy > 0.0f ? a.occupancy : 1.0f);
    p = pf_fixed<6, 100, 2>(p, a.tempFactor);
    for (int i = 0; i < 10; ++i) *p++ = ' ';
    *p++ = ' ';
    *p++ = (a.atom.data[0] == '\0') ? ' ' : a.atom.data[0];
    *p++ = ' '; *p++ = ' '; *p++ = '\n';
    // total: 81 bytes
}
