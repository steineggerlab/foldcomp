/**
 * File: atom_coordinate.cpp
 * Project: foldcomp
 * Created: 2021-01-18 12:53:34
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     The data type to handle atom coordinate comes here.
 * ---
 * Last Modified: Fri Mar 03 2023
 * Modified By: Hyunbin Kim
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "atom_coordinate.h"
#include "amino_acid.h"
#include "foldcomp.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
#include "gemmi/model.hpp"
#define GEMMI_WRITE_IMPLEMENTATION
#include "gemmi/to_mmcif.hpp"
#include "gemmi/to_cif.hpp"
#include <sstream>
#endif

namespace {

constexpr float MAX_PEPTIDE_BOND_DISTANCE = 2.5f;
constexpr float MAX_PEPTIDE_BOND_DISTANCE_SQUARED =
    MAX_PEPTIDE_BOND_DISTANCE * MAX_PEPTIDE_BOND_DISTANCE;

void appendSpaces(std::string& output, size_t count) {
    output.append(count, ' ');
}

void appendLeftAligned(std::string& output, const std::string& value, size_t width) {
    size_t copied = std::min(width, value.size());
    output.append(value.data(), copied);
    if (copied < width) {
        appendSpaces(output, width - copied);
    }
}

void appendRightAligned(std::string& output, const std::string& value, size_t width) {
    size_t copied = std::min(width, value.size());
    if (copied < width) {
        appendSpaces(output, width - copied);
    }
    output.append(value.data() + (value.size() - copied), copied);
}

void appendRightAlignedInt(std::string& output, int value, size_t width) {
    bool negative = value < 0;
    uint32_t magnitude;
    if (negative) {
        magnitude = static_cast<uint32_t>(-(static_cast<int64_t>(value)));
    } else {
        magnitude = static_cast<uint32_t>(value);
    }
    char digits[16];
    int len = 0;
    do {
        digits[len++] = static_cast<char>('0' + (magnitude % 10));
        magnitude /= 10;
    } while (magnitude != 0);
    size_t totalLen = static_cast<size_t>(len + (negative ? 1 : 0));
    if (totalLen < width) {
        appendSpaces(output, width - totalLen);
    }
    if (negative) {
        output.push_back('-');
    }
    while (len-- > 0) {
        output.push_back(digits[len]);
    }
}

template <size_t Width, uint32_t Scale, size_t Precision>
void appendRightAlignedFixed(std::string& output, float value) {
    bool negative = value < 0.0f;
    float adjusted = value + (negative ? -(0.5f / static_cast<float>(Scale))
                                       :  (0.5f / static_cast<float>(Scale)));
    int64_t scaled = static_cast<int64_t>(adjusted * static_cast<float>(Scale));
    if (scaled < 0) {
        scaled = -scaled;
    }
    uint32_t fraction = static_cast<uint32_t>(scaled % Scale);
    uint32_t integer = static_cast<uint32_t>(scaled / Scale);

    char intDigits[16];
    int intLen = 0;
    do {
        intDigits[intLen++] = static_cast<char>('0' + (integer % 10));
        integer /= 10;
    } while (integer != 0);

    size_t totalLen = static_cast<size_t>(intLen) + 1 + Precision + (negative ? 1 : 0);
    if (totalLen < Width) {
        appendSpaces(output, Width - totalLen);
    }
    if (negative) {
        output.push_back('-');
    }
    while (intLen-- > 0) {
        output.push_back(intDigits[intLen]);
    }
    output.push_back('.');

    char fracDigits[Precision];
    for (size_t i = 0; i < Precision; i++) {
        fracDigits[Precision - 1 - i] = static_cast<char>('0' + (fraction % 10));
        fraction /= 10;
    }
    output.append(fracDigits, Precision);
}

void appendTitleLines(std::string& output, const std::string& title) {
    if (title.empty()) {
        return;
    }
    size_t offset = 0;
    int continuation = 1;
    while (offset < title.size()) {
        size_t chunk = std::min<size_t>(70, title.size() - offset);
        if (continuation == 1) {
            output.append("TITLE     ", 10);
        } else {
            output.append("TITLE  ", 7);
            appendRightAlignedInt(output, continuation, 3);
        }
        output.append(title.data() + offset, chunk);
        output.push_back('\n');
        offset += chunk;
        continuation++;
    }
}

bool residueRequiresRawEncoding(const tcb::span<const AtomCoordinate>& residueAtoms) {
    for (const auto& atom : residueAtoms) {
        if (atom.altloc != ' ' && atom.altloc != '\0') {
            return true;
        }
        if (atom.insertion_code != ' ' && atom.insertion_code != '\0') {
            return true;
        }
    }
    return false;
}

struct ResidueSummary {
    size_t start;
    size_t end;
    int residueIndex;
    const AtomCoordinate* n;
    const AtomCoordinate* ca;
    const AtomCoordinate* c;
    bool requiresRawEncoding;
};

bool hasAbnormalPeptideGap(const ResidueSummary& previousResidue, const ResidueSummary& currentResidue) {
    if (previousResidue.c == nullptr || currentResidue.n == nullptr) {
        return false;
    }
    float dx = previousResidue.c->coordinate.x - currentResidue.n->coordinate.x;
    float dy = previousResidue.c->coordinate.y - currentResidue.n->coordinate.y;
    float dz = previousResidue.c->coordinate.z - currentResidue.n->coordinate.z;
    float distanceSquared = dx * dx + dy * dy + dz * dz;
    return distanceSquared > MAX_PEPTIDE_BOND_DISTANCE_SQUARED;
}

std::vector<ResidueSummary> buildResidueSummaries(
    const tcb::span<const AtomCoordinate>& atoms,
    size_t offset = 0
) {
    std::vector<ResidueSummary> summaries;
    if (atoms.empty()) {
        return summaries;
    }

    summaries.reserve(atoms.size() / 4 + 1);
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) || startsNewResidue(atoms[i], atoms[i - 1]);
        if (!endOfResidue) {
            continue;
        }

        const AtomCoordinate* n = nullptr;
        const AtomCoordinate* ca = nullptr;
        const AtomCoordinate* c = nullptr;
        for (size_t j = residueStart; j < i; j++) {
            const AtomCoordinate& atom = atoms[j];
            if (atom.atom == "N" && n == nullptr) {
                n = &atom;
            } else if (atom.atom == "CA" && ca == nullptr) {
                ca = &atom;
            } else if (atom.atom == "C" && c == nullptr) {
                c = &atom;
            }
        }
        summaries.push_back({
            residueStart + offset, i + offset, atoms[residueStart].residue_index,
            n, ca, c,
            residueRequiresRawEncoding(tcb::span<const AtomCoordinate>(atoms.data() + residueStart, i - residueStart))
        });
        residueStart = i;
    }

    return summaries;
}

#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
std::string sanitizeMMCIFDataName(const std::string& title) {
    if (title.empty()) {
        return "foldcomp";
    }
    std::string output;
    output.reserve(title.size());
    for (char ch : title) {
        unsigned char uch = static_cast<unsigned char>(ch);
        if (std::isalnum(uch) || ch == '_') {
            output.push_back(ch);
        } else {
            output.push_back('_');
        }
    }
    if (output.empty() || output[0] == '_' || output[0] == '#') {
        output.insert(output.begin(), 'd');
    }
    return output;
}

gemmi::Element inferGemmiElement(const std::string& atomName) {
    const char* ptr = atomName.c_str();
    while (*ptr != '\0' && !std::isalpha(static_cast<unsigned char>(*ptr))) {
        ++ptr;
    }
    return gemmi::find_element(ptr);
}
#endif

}

/**
 * @brief Construct a new Atom Coordinate:: Atom Coordinate object
 *
 * @param a A string for atom name
 * @param r A string for residue name
 * @param ai An integer for atom index
 * @param ri An integer for residue index
 * @param x A float for x coordinate
 * @param y A float for y coordinate
 * @param z A float for z coordinate
 */
AtomCoordinate::AtomCoordinate(
    const std::string& a, const std::string& r, const std::string& c,
    int ai, int ri, float x, float y, float z,
    float occupancy, float tempFactor, int model, char insertionCode, char altloc
): atom(a), residue(r), chain(c), atom_index(ai), residue_index(ri), occupancy(occupancy), tempFactor(tempFactor), model(model), insertion_code(insertionCode), altloc(altloc) {
    this->coordinate = {x, y, z};
}

/**
 * @brief Construct a new Atom Coordinate:: Atom Coordinate object
 *
 * @param a A string for atom name
 * @param r A string for residue name
 * @param ai An integer for atom index
 * @param ri An integer for residue index
 * @param coord A float vector for x,y,z coordinates.
 */
AtomCoordinate::AtomCoordinate(
    const std::string& a, const std::string& r, const std::string& c,
    int ai, int ri, float3d coord,
    float occupancy, float tempFactor, int model, char insertionCode, char altloc
): atom(a), residue(r), chain(c), atom_index(ai), residue_index(ri), coordinate(coord), occupancy(occupancy), tempFactor(tempFactor), model(model), insertion_code(insertionCode), altloc(altloc) {
}

bool startsNewResidue(const AtomCoordinate& current, const AtomCoordinate& previous) {
    if (current.chain != previous.chain) {
        return true;
    }
    if (current.model != previous.model) {
        return true;
    }
    if (current.residue_index != previous.residue_index) {
        return true;
    }
    if (current.insertion_code != previous.insertion_code) {
        return true;
    }
    if (current.residue != previous.residue) {
        return true;
    }
    if (current.atom == "N" && previous.atom != "N") {
        // Reused residue indices can still denote a new residue in author numbering.
        return true;
    }
    return false;
}

bool AtomCoordinate::operator==(const AtomCoordinate& other) const {
    return (
        (this->atom == other.atom) &&
        (this->atom_index == other.atom_index) &&
        (this->residue == other.residue) &&
        (this->residue_index == other.residue_index) &&
        (this->chain == other.chain) &&
        (this->model == other.model) &&
        (this->insertion_code == other.insertion_code) &&
        (this->altloc == other.altloc) &&
        (this->coordinate.x == other.coordinate.x) &&
        (this->coordinate.y == other.coordinate.y) &&
        (this->coordinate.z == other.coordinate.z)
    );
}
bool AtomCoordinate::operator!=(const AtomCoordinate& other) const {
    return !(*this == other);
}

bool AtomCoordinate::isBackbone() const {
    return ((this->atom == "N") ||(this->atom == "CA") ||(this->atom == "C"));
}

void AtomCoordinate::print(int option) const {
    std::cout << "Atom: " << this->atom << std::endl;
    if (option != 0) {
        std::cout << "Residue: " << this->residue << std::endl;
        std::cout << "Chain: " << this->chain << std::endl;
        std::cout << "Atom Index: " << this->atom_index << std::endl;
        std::cout << "Residue Index: " << this->residue_index << std::endl;
        if (option == 2) {
            std::cout << "Coordinate: ";
            std::cout << this->coordinate.x << " ";
            std::cout << this->coordinate.y << " ";
            std::cout << this->coordinate.z << " ";
            std::cout << std::endl;
        }
    }
}

/**
 * @brief Extracts coordinates from AtomCoordinate vector
 *
 * @param atoms A vector of AtomCoordinate
 * @return std::vector< std::vector<float> >
 */
std::vector<float3d> extractCoordinates(
    const std::vector<AtomCoordinate>& atoms
) {
    std::vector<float3d> output(atoms.size());
    for (size_t i = 0; i < atoms.size(); i++) {
        output[i] = atoms[i].coordinate;
    }
    return output;
}

std::vector<AtomCoordinate> extractChain(
    std::vector<AtomCoordinate>& atoms, std::string chain
) {
    std::vector<AtomCoordinate> output;
    int total = atoms.size();
    for (int i = 0; i < total; i++) {
        const AtomCoordinate& curr_atm = atoms[i];
        if (i < (total-1)) {
            const AtomCoordinate& next_atm = atoms[i + 1];
            if (next_atm.atom == curr_atm.atom) {
                continue;
            }
        }
        if (curr_atm.chain == chain) {
            output.push_back(curr_atm);
        }
    }
    return output;
}

void printAtomCoordinateVector(std::vector<AtomCoordinate>& atoms, int option) {
    for (const AtomCoordinate& curr_atm : atoms) {
        curr_atm.print(option);
    }
}

std::vector<AtomCoordinate> filterBackbone(const tcb::span<AtomCoordinate>& atoms) {
    std::vector<AtomCoordinate> output;
    std::vector<std::pair<size_t, size_t>> residueRanges = splitResidueRanges(atoms);
    output.reserve(residueRanges.size() * 3);
    for (const auto& residueRange : residueRanges) {
        const AtomCoordinate* n = nullptr;
        const AtomCoordinate* ca = nullptr;
        const AtomCoordinate* c = nullptr;
        for (size_t i = residueRange.first; i < residueRange.second; i++) {
            const auto& atom = atoms[i];
            if (atom.atom == "N" && n == nullptr) {
                n = &atom;
            } else if (atom.atom == "CA" && ca == nullptr) {
                ca = &atom;
            } else if (atom.atom == "C" && c == nullptr) {
                c = &atom;
            }
        }
        if (n != nullptr && ca != nullptr && c != nullptr) {
            output.emplace_back(*n);
            output.emplace_back(*ca);
            output.emplace_back(*c);
        }
    }
    return output;
}

std::vector<AtomCoordinate> weightedAverage(
    const std::vector<AtomCoordinate>& origAtoms, const std::vector<AtomCoordinate>& revAtoms
) {
    std::vector<AtomCoordinate> output;
    output.reserve(origAtoms.size());
    int total = origAtoms.size();
    for (int i = 0; i < total; i++) {
        const AtomCoordinate& curr_atm = origAtoms[i];
        const AtomCoordinate&  rev_atm = revAtoms[i];
        output.emplace_back(
            curr_atm.atom, curr_atm.residue, curr_atm.chain,
            curr_atm.atom_index, curr_atm.residue_index,
            ((curr_atm.coordinate.x * (float)(total - i)) + (rev_atm.coordinate.x * (float)i)) / (float)total,
            ((curr_atm.coordinate.y * (float)(total - i)) + (rev_atm.coordinate.y * (float)i)) / (float)total,
            ((curr_atm.coordinate.z * (float)(total - i)) + (rev_atm.coordinate.z * (float)i)) / (float)total
        );
    }
    return output;
}

void writeAtomCoordinatesToPDB(
    const std::vector<AtomCoordinate>& atoms, const std::string& title, std::string& output,
    bool appendOutput, bool emitFinalTer
) {
    if (!appendOutput) {
        output.clear();
    }
    output.reserve(output.size() + title.size() + atoms.size() * 96);
    if (!title.empty()) {
        appendTitleLines(output, title);
    }

    int total = static_cast<int>(atoms.size());
    for (int i = 0; i < total; i++) {
        const AtomCoordinate& atom = atoms[i];
        output.append("ATOM  ", 6);
        appendRightAlignedInt(output, atom.atom_index, 5);
        output.push_back(' ');
        if (atom.atom.size() == 4) {
            appendLeftAligned(output, atom.atom, 4);
        } else {
            output.push_back(' ');
            appendLeftAligned(output, atom.atom, 3);
        }
        output.push_back(atom.altloc == '\0' ? ' ' : atom.altloc);
        appendRightAligned(output, atom.residue, 3);
        output.push_back(' ');
        char chainId = atom.chain.empty() ? ' ' : atom.chain[0];
        output.push_back(chainId);
        appendRightAlignedInt(output, atom.residue_index, 4);
        output.push_back(atom.insertion_code == '\0' ? ' ' : atom.insertion_code);
        output.append("   ", 3);
        appendRightAlignedFixed<8, 1000, 3>(output, atom.coordinate.x);
        appendRightAlignedFixed<8, 1000, 3>(output, atom.coordinate.y);
        appendRightAlignedFixed<8, 1000, 3>(output, atom.coordinate.z);
        appendRightAlignedFixed<6, 100, 2>(output, atom.occupancy > 0.0f ? atom.occupancy : 1.0f);
        appendRightAlignedFixed<6, 100, 2>(output, atom.tempFactor);
        output.append("          ", 10);
        output.push_back(' ');
        output.push_back(atom.atom.empty() ? ' ' : atom.atom[0]);
        output.append("  \n", 3);
        bool needsTer = false;
        if (i == (total - 1)) {
            needsTer = emitFinalTer;
        } else {
            const AtomCoordinate& nextAtom = atoms[i + 1];
            needsTer = (nextAtom.model != atom.model) || (nextAtom.chain != atom.chain);
        }
        if (needsTer) {
            output.append("TER   ", 6);
            appendRightAlignedInt(output, atom.atom_index + 1, 5);
            output.append("      ", 6);
            appendRightAligned(output, atom.residue, 3);
            output.push_back(' ');
            output.push_back(chainId);
            appendRightAlignedInt(output, atom.residue_index, 4);
            output.push_back(atom.insertion_code == '\0' ? ' ' : atom.insertion_code);
            output.push_back('\n');
        }
    }
}

int writeAtomCoordinatesToPDBFile(
    const std::vector<AtomCoordinate>& atoms, const std::string& title, const std::string& pdb_path
) {
    std::string output;
    writeAtomCoordinatesToPDB(atoms, title, output);
    FILE* pdb_file = fopen(pdb_path.c_str(), "wb");
    if (pdb_file == nullptr) {
        return 1;
    }
    size_t written = fwrite(output.data(), 1, output.size(), pdb_file);
    fclose(pdb_file);
    return written == output.size() ? 0 : 1;
}

#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
bool writeAtomCoordinatesToMMCIF(
    std::vector<AtomCoordinate>& atoms, const std::string& title, std::string& output
) {
    gemmi::Structure st;
    st.name = sanitizeMMCIFDataName(title);
    if (!title.empty()) {
        st.info["_entry.id"] = st.name;
        st.info["_struct.title"] = title;
    }

    std::vector<std::pair<size_t, size_t>> residueRanges = splitResidueRanges(tcb::span<AtomCoordinate>(atoms.data(), atoms.size()));
    int currentModel = 0;
    std::string currentChain;
    int currentLabelSeq = 0;
    gemmi::Model* model = nullptr;
    gemmi::Chain* chain = nullptr;

    for (const auto& residueRange : residueRanges) {
        const AtomCoordinate& firstAtom = atoms[residueRange.first];
        if (model == nullptr || firstAtom.model != currentModel) {
            st.models.emplace_back(std::to_string(firstAtom.model));
            model = &st.models.back();
            chain = nullptr;
            currentModel = firstAtom.model;
            currentChain.clear();
            currentLabelSeq = 0;
        }
        if (chain == nullptr || firstAtom.chain != currentChain) {
            model->chains.emplace_back(firstAtom.chain);
            chain = &model->chains.back();
            currentChain = firstAtom.chain;
            currentLabelSeq = 0;
        }
        currentLabelSeq++;

        gemmi::Residue residue;
        residue.name = firstAtom.residue;
        residue.seqid = gemmi::SeqId(firstAtom.residue_index, firstAtom.insertion_code == '\0' ? ' ' : firstAtom.insertion_code);
        residue.subchain = currentChain;
        residue.entity_id = "1";
        residue.label_seq = currentLabelSeq;
        residue.entity_type = gemmi::EntityType::Polymer;
        residue.het_flag = 'A';

        for (size_t atomIndex = residueRange.first; atomIndex < residueRange.second; atomIndex++) {
            const AtomCoordinate& atomCoordinate = atoms[atomIndex];
            gemmi::Atom atom;
            atom.name = atomCoordinate.atom;
            atom.serial = atomCoordinate.atom_index;
            atom.pos = gemmi::Position(atomCoordinate.coordinate.x, atomCoordinate.coordinate.y, atomCoordinate.coordinate.z);
            atom.occ = atomCoordinate.occupancy > 0.0f ? atomCoordinate.occupancy : 1.0f;
            atom.b_iso = atomCoordinate.tempFactor;
            atom.element = inferGemmiElement(atomCoordinate.atom);
            atom.altloc = (atomCoordinate.altloc == '\0' || atomCoordinate.altloc == ' ')
                              ? '\0'
                              : atomCoordinate.altloc;
            residue.atoms.push_back(std::move(atom));
        }
        chain->residues.push_back(std::move(residue));
    }

    gemmi::MmcifOutputGroups groups(false);
    groups.atoms = true;
    groups.block_name = true;
    groups.entry = true;
    groups.title_keywords = !title.empty();
    groups.struct_asym = true;
    groups.group_pdb = true;

    gemmi::cif::Document doc = gemmi::make_mmcif_document(st, groups);
    std::ostringstream stream;
    gemmi::cif::write_cif_to_stream(stream, doc, gemmi::cif::Style::Pdbx);
    output = stream.str();
    return true;
}
#endif

std::vector< std::vector<AtomCoordinate> > splitAtomByResidue(
    const tcb::span<AtomCoordinate>& atomCoordinates
) {
    std::vector< std::vector<AtomCoordinate> > output;
    std::vector<std::pair<size_t, size_t>> ranges = splitResidueRanges(atomCoordinates);
    output.reserve(ranges.size());
    for (const auto& range : ranges) {
        output.emplace_back(atomCoordinates.begin() + range.first, atomCoordinates.begin() + range.second);
    }
    return output;
}

std::vector<std::pair<size_t, size_t>> splitResidueRanges(
    const tcb::span<AtomCoordinate>& atomCoordinates
) {
    std::vector<std::pair<size_t, size_t>> output;
    if (atomCoordinates.empty()) {
        return output;
    }
    output.reserve(atomCoordinates.size() / 4 + 1);
    size_t start = 0;
    for (size_t i = 1; i < atomCoordinates.size(); i++) {
        if (startsNewResidue(atomCoordinates[i], atomCoordinates[i - 1])) {
            output.push_back({start, i});
            start = i;
        }
    }
    output.push_back({start, atomCoordinates.size()});
    return output;
}

std::vector<std::string> getResidueNameVector(
    const tcb::span<AtomCoordinate>& atomCoordinates
) {
    std::vector<std::string> output;
    std::vector<std::pair<size_t, size_t>> residueRanges = splitResidueRanges(atomCoordinates);
    output.reserve(residueRanges.size());
    for (const auto& residueRange : residueRanges) {
        output.push_back(atomCoordinates[residueRange.first].residue);
    }
    return output;
}

AtomCoordinate findFirstAtom(const std::vector<AtomCoordinate>& atoms, std::string atom_name) {
    return findFirstAtom(tcb::span<const AtomCoordinate>(atoms.data(), atoms.size()), std::move(atom_name));
}

AtomCoordinate findFirstAtom(const tcb::span<const AtomCoordinate>& atoms, std::string atom_name) {
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.atom == atom_name) {
            return curr_atm;
        }
    }
    return AtomCoordinate();
}

void setAtomIndexSequentially(std::vector<AtomCoordinate>& atoms, int start) {
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms[i].atom_index = start + i;
    }
}

void removeAlternativePosition(std::vector<AtomCoordinate>& atoms) {
    if (atoms.empty()) {
        return;
    }
    size_t writeIndex = 1;
    for (size_t readIndex = 1; readIndex < atoms.size(); readIndex++) {
        const AtomCoordinate& current = atoms[readIndex];
        const AtomCoordinate& previous = atoms[writeIndex - 1];
        if (current.atom == previous.atom &&
            current.residue == previous.residue &&
            current.residue_index == previous.residue_index &&
            current.chain == previous.chain &&
            current.model == previous.model &&
            current.insertion_code == previous.insertion_code) {
            continue;
        }
        if (writeIndex != readIndex) {
            atoms[writeIndex] = std::move(atoms[readIndex]);
        }
        writeIndex++;
    }
    atoms.resize(writeIndex);
}

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    std::vector<AtomCoordinate>& atoms, int residue_index
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index == residue_index) {
            output.emplace_back(curr_atm);
        }
    }
    return output;
}

std::vector<AtomCoordinate> getAtomsWithResidueIndiceRange(
    std::vector<AtomCoordinate>& atoms, int start, int end
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index >= start && curr_atm.residue_index < end) {
            output.emplace_back(curr_atm);
        }
    }
    return output;
}

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, int residue_index,
    std::vector<std::string> atomNames
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index == residue_index) {
            for (const std::string& atom_name : atomNames) {
                if (curr_atm.atom == atom_name) {
                    output.emplace_back(curr_atm);
                }
            }
        }
    }
    return output;
}

std::vector< std::vector<AtomCoordinate> > getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, std::vector<int> residue_index,
    std::vector<std::string> atomNames
) {
    std::vector<std::vector<AtomCoordinate>> output;
    for (int curr_index : residue_index) {
        output.emplace_back(getAtomsWithResidueIndex(atoms, curr_index, atomNames));
    }
    return output;
}

float RMSD(std::vector<AtomCoordinate>& atoms1, std::vector<AtomCoordinate>& atoms2) {
    // RMSD: Root Mean Square Deviation
    float sum = 0;
    // Sum of square of distance
    for (size_t i = 0; i < atoms1.size(); i++) {
        sum += pow(atoms1[i].coordinate.x - atoms2[i].coordinate.x, 2);
        sum += pow(atoms1[i].coordinate.y - atoms2[i].coordinate.y, 2);
        sum += pow(atoms1[i].coordinate.z - atoms2[i].coordinate.z, 2);
    }
    return sqrt(sum / atoms1.size());
}

std::vector<AtomCoordinate> _subsetAtomVectorWithIndices(
    std::vector<AtomCoordinate>& atoms,
    std::pair<size_t, size_t>& indices
) {
    std::vector<AtomCoordinate> output;
    for (size_t i = indices.first; i < indices.second; i++) {
        output.push_back(atoms[i]);
    }
    return output;
}

void _splitAtomVectorWithIndices(
    std::vector<AtomCoordinate>& atoms,
    std::vector< std::pair<size_t, size_t> >& indices,
    std::vector< std::vector<AtomCoordinate> >& output
) {
    if (indices.size() == 0) {
        output.push_back(atoms);
    } else {
        for (size_t i = 0; i < indices.size(); i++) {
            output.push_back(_subsetAtomVectorWithIndices(atoms, indices[i]));
        }
    }
}

/**
 * @brief Identify discontinuous regions in atom coordinate vector and return
 *        vector of indices of the start and end of each region.
 *        start: inclusive, end: exclusive [start, end)
 * @param atoms
 * @param mode
 * @return std::vector<std::pair<size_t, size_t>>
 */
std::vector< std::pair<size_t, size_t> > identifyChains(const std::vector<AtomCoordinate>& atoms) {
    std::vector< std::pair<size_t, size_t> > output;
    if (atoms.empty()) {
        return output;
    }
    size_t start = 0;
    for (size_t i = 1; i < atoms.size(); i++) {
        if (atoms[i].model != atoms[i - 1].model || atoms[i].chain != atoms[i - 1].chain) {
            output.emplace_back(start, i);
            start = i;
        }
    }
    if (start < atoms.size()) {
        output.emplace_back(start, atoms.size());
    }
    return output;
}

/**
 * @brief Identify discontinuous residue indices in atom coordinate vector and return
 *        coordinates that have the same chain
 * @param atoms
 * @return std::vector< std::pair<size_t, size_t> >
 */
std::vector<std::pair<size_t, size_t>> identifyDiscontinousResInd(
    const std::vector<AtomCoordinate>& atoms,
    size_t chain_start,
    size_t chain_end
) {
    std::vector<std::pair<size_t, size_t>> output;
    if (chain_start >= chain_end || chain_end > atoms.size()) {
        return output;
    }
    tcb::span<const AtomCoordinate> chainSpan(atoms.data() + chain_start, chain_end - chain_start);
    std::vector<ResidueSummary> summaries = buildResidueSummaries(chainSpan, chain_start);
    if (summaries.empty()) {
        return output;
    }
    size_t start = chain_start;
    for (size_t i = 1; i < summaries.size(); i++) {
        if (summaries[i].residueIndex != summaries[i - 1].residueIndex + 1 ||
            hasAbnormalPeptideGap(summaries[i - 1], summaries[i])) {
            output.emplace_back(start, summaries[i].start);
            start = summaries[i].start;
        }
    }
    output.emplace_back(start, chain_end);
    return output;
}

std::vector<std::pair<size_t, size_t>> identifyCompleteBackboneRegions(
    const tcb::span<AtomCoordinate>& atoms
) {
    std::vector<std::pair<size_t, size_t>> output;
    if (atoms.empty()) {
        return output;
    }

    std::vector<ResidueSummary> residues = buildResidueSummaries(
        tcb::span<const AtomCoordinate>(atoms.data(), atoms.size())
    );

    size_t i = 0;
    while (i < residues.size()) {
        if (!(residues[i].n != nullptr && residues[i].ca != nullptr && residues[i].c != nullptr &&
              !residues[i].requiresRawEncoding)) {
            i++;
            continue;
        }
        size_t runStartIndex = i;
        size_t start = residues[i].start;
        size_t end = residues[i].end;
        i++;
        while (i < residues.size() &&
               residues[i].n != nullptr && residues[i].ca != nullptr && residues[i].c != nullptr &&
               !residues[i].requiresRawEncoding) {
            end = residues[i].end;
            i++;
        }
        if (i - runStartIndex >= 2) {
            output.emplace_back(start, end);
        }
    }
    return output;
}

std::vector<BackboneRegion> identifyBackboneRegions(const tcb::span<AtomCoordinate>& atoms) {
    std::vector<BackboneRegion> output;
    if (atoms.empty()) {
        return output;
    }

    std::vector<ResidueSummary> residues = buildResidueSummaries(
        tcb::span<const AtomCoordinate>(atoms.data(), atoms.size())
    );

    size_t i = 0;
    while (i < residues.size()) {
        bool completeBackbone = residues[i].n != nullptr && residues[i].ca != nullptr &&
                                residues[i].c != nullptr && !residues[i].requiresRawEncoding;
        if (!completeBackbone) {
            size_t start = residues[i].start;
            size_t end = residues[i].end;
            i++;
            while (i < residues.size()) {
                bool nextCompleteBackbone = residues[i].n != nullptr && residues[i].ca != nullptr &&
                                            residues[i].c != nullptr && !residues[i].requiresRawEncoding;
                if (nextCompleteBackbone) {
                    break;
                }
                end = residues[i].end;
                i++;
            }
            output.push_back({start, end, false});
            continue;
        }

        size_t runStartIndex = i;
        size_t start = residues[i].start;
        size_t end = residues[i].end;
        i++;
        while (i < residues.size() &&
               residues[i].n != nullptr && residues[i].ca != nullptr && residues[i].c != nullptr &&
               !residues[i].requiresRawEncoding) {
            end = residues[i].end;
            i++;
        }
        const size_t runLength = i - runStartIndex;
        if (runLength >= 2) {
            output.push_back({start, end, true});
        } else {
            output.push_back({start, end, false});
        }
    }

    return output;
}

namespace {

template <typename T>
void appendBinaryValue(std::string& output, const T& value) {
    size_t offset = output.size();
    output.resize(offset + sizeof(T));
    memcpy(&output[offset], &value, sizeof(T));
}

template <typename T>
bool readBinaryValue(const char* data, size_t size, size_t& offset, T& out) {
    if (offset + sizeof(T) > size) {
        return false;
    }
    memcpy(&out, data + offset, sizeof(T));
    offset += sizeof(T);
    return true;
}

} // namespace

bool serializeAtomCoordinates(
    const tcb::span<const AtomCoordinate>& atoms,
    std::string& output
) {
    output.clear();
    output.reserve(sizeof(uint32_t) + atoms.size() * 40);
    uint32_t atomCount = static_cast<uint32_t>(atoms.size());
    appendBinaryValue(output, atomCount);
    for (const auto& atom : atoms) {
        if (atom.atom.size() > 255 || atom.residue.size() > 255) {
            return false;
        }
        uint8_t atomNameLen = static_cast<uint8_t>(atom.atom.size());
        uint8_t residueNameLen = static_cast<uint8_t>(atom.residue.size());
        appendBinaryValue(output, atomNameLen);
        appendBinaryValue(output, residueNameLen);
        appendBinaryValue(output, atom.atom_index);
        appendBinaryValue(output, atom.residue_index);
        appendBinaryValue(output, atom.coordinate.x);
        appendBinaryValue(output, atom.coordinate.y);
        appendBinaryValue(output, atom.coordinate.z);
        appendBinaryValue(output, atom.occupancy);
        appendBinaryValue(output, atom.tempFactor);
        appendBinaryValue(output, atom.insertion_code);
        appendBinaryValue(output, atom.altloc);
        output.append(atom.atom);
        output.append(atom.residue);
    }
    return true;
}

bool serializeAtomCoordinates(
    const std::vector<AtomCoordinate>& atoms,
    std::string& output
) {
    return serializeAtomCoordinates(tcb::span<const AtomCoordinate>(atoms.data(), atoms.size()), output);
}

bool deserializeAtomCoordinates(
    const char* data,
    size_t size,
    std::vector<AtomCoordinate>& atoms
) {
    auto parse = [&](bool hasIdentityFields) -> bool {
        atoms.clear();
        size_t offset = 0;
        uint32_t atomCount = 0;
        if (!readBinaryValue(data, size, offset, atomCount)) {
            return false;
        }
        atoms.reserve(atomCount);
        for (uint32_t i = 0; i < atomCount; i++) {
            uint8_t atomNameLen = 0;
            uint8_t residueNameLen = 0;
            AtomCoordinate atom;
            if (!readBinaryValue(data, size, offset, atomNameLen) ||
                !readBinaryValue(data, size, offset, residueNameLen) ||
                !readBinaryValue(data, size, offset, atom.atom_index) ||
                !readBinaryValue(data, size, offset, atom.residue_index) ||
                !readBinaryValue(data, size, offset, atom.coordinate.x) ||
                !readBinaryValue(data, size, offset, atom.coordinate.y) ||
                !readBinaryValue(data, size, offset, atom.coordinate.z) ||
                !readBinaryValue(data, size, offset, atom.occupancy) ||
                !readBinaryValue(data, size, offset, atom.tempFactor)) {
                return false;
            }
            if (hasIdentityFields) {
                if (!readBinaryValue(data, size, offset, atom.insertion_code) ||
                    !readBinaryValue(data, size, offset, atom.altloc)) {
                    return false;
                }
            } else {
                atom.insertion_code = ' ';
                atom.altloc = ' ';
            }
            size_t required = static_cast<size_t>(atomNameLen) + static_cast<size_t>(residueNameLen);
            if (offset + required > size) {
                return false;
            }
            atom.atom.assign(data + offset, data + offset + atomNameLen);
            offset += atomNameLen;
            atom.residue.assign(data + offset, data + offset + residueNameLen);
            offset += residueNameLen;
            atoms.push_back(std::move(atom));
        }
        return offset == size;
    };
    if (parse(true)) {
        return true;
    }
    return parse(false);
}
