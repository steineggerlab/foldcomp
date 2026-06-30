#include "amino_acid.h"
#include "atom_coordinate.h"
#include "database_reader.h"
#include "foldcomp.h"
#include "structure_codec.h"
#include "structure_reader.h"
#include "torsion_angle.h"
#include "utility.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct ResidueStats {
    std::set<std::string> atoms;
    std::string residue;
};

struct Summary {
    size_t atomCount = 0;
    size_t residueCount = 0;
    size_t completeBackboneResidueCount = 0;
    size_t caOnlyResidueCount = 0;
};

using ResidueKey = std::tuple<int, std::string, int, int>;
using RmsdResidueKey = std::tuple<int, std::string, int, int>;
namespace fs = std::filesystem;

struct ParsedFczPayload {
    CompressedFileHeader header{};
    std::vector<int> anchorIndices;
    std::string title;
    std::vector<std::array<float, 3>> prevCoords;
    std::vector<std::array<float, 3>> innerAnchorCoords;
    std::vector<std::array<float, 3>> lastCoords;
    char hasOXT = 0;
    std::array<float, 3> oxtCoords{};
    std::vector<BackboneChain> compressedBackbone;
    std::vector<unsigned char> sideChainAngles;
    float tempMin = 0.0f;
    float tempContF = 0.0f;
    std::vector<unsigned char> tempFactors;
};

struct HeaderFieldDiff {
    std::string name;
    std::string leftValue;
    std::string rightValue;
};

struct TraceFragmentSelection {
    int globalFragmentIndex = -1;
    size_t chainRegionIndex = 0;
    size_t discontinuityIndex = 0;
    size_t regionIndex = 0;
    bool useFoldcompEncoding = false;
    std::vector<AtomCoordinate> regionAtoms;
};

Summary summarize(const std::vector<AtomCoordinate>& atoms);
std::map<ResidueKey, ResidueStats> groupResidues(const std::vector<AtomCoordinate>& atoms);
std::map<RmsdResidueKey, std::vector<AtomCoordinate>> groupResiduesForRmsd(const std::vector<AtomCoordinate>& atoms);
bool calculateAlignedRmsd(
    const std::vector<AtomCoordinate>& atoms1,
    const std::vector<AtomCoordinate>& atoms2,
    float& backboneRmsd,
    float& allAtomRmsd
);
bool parseFczPayload(const std::string& payload, ParsedFczPayload& parsed);
bool firstHeaderFieldDiff(const ParsedFczPayload& left, const ParsedFczPayload& right, HeaderFieldDiff& diff);
bool findTraceFragment(const std::vector<AtomCoordinate>& atoms, int selectedFragmentIndex, TraceFragmentSelection& selection);

enum class RawReason {
    NonProtein = 0,
    CaOnlyResidue,
    IncompleteBackboneResidue,
    IsolatedCompleteResidue,
    UnknownResidue,
    PartialCanonicalAtomSet,
    ExtraNonCanonicalAtoms,
};

struct FallbackStats {
    size_t fileCount = 0;
    size_t proteinFileCount = 0;
    size_t rawFileCount = 0;
    size_t residueCount = 0;
    std::map<RawReason, size_t> rawResiduesByReason;
    std::map<RawReason, size_t> rawFilesByReason;
};

const char* rawReasonName(RawReason reason) {
    switch (reason) {
        case RawReason::NonProtein: return "non_protein";
        case RawReason::CaOnlyResidue: return "ca_only_residue";
        case RawReason::IncompleteBackboneResidue: return "incomplete_backbone_residue";
        case RawReason::IsolatedCompleteResidue: return "isolated_complete_residue";
        case RawReason::UnknownResidue: return "unknown_residue";
        case RawReason::PartialCanonicalAtomSet: return "partial_canonical_atom_set";
        case RawReason::ExtraNonCanonicalAtoms: return "extra_noncanonical_atoms";
    }
    return "unknown";
}

bool residueNeedsRawFallbackForTrace(const tcb::span<const AtomCoordinate>& residueAtoms) {
    if (residueAtoms.empty()) {
        return false;
    }
    for (const auto& atom : residueAtoms) {
        if (atom.altloc != ' ' && atom.altloc != '\0') {
            return true;
        }
        if (atom.insertion_code != ' ' && atom.insertion_code != '\0') {
            return true;
        }
    }
    auto aaIt = Foldcomp::AAS.find(residueAtoms[0].residue);
    if (aaIt == Foldcomp::AAS.end()) {
        return true;
    }
    std::set<std::string> sourceAtoms;
    for (const auto& atom : residueAtoms) {
        sourceAtoms.insert(atom.atom);
    }
    std::set<std::string> canonicalAtoms(aaIt->second.atoms.begin(), aaIt->second.atoms.end());
    std::set<std::string> canonicalAtomsWithOxt = canonicalAtoms;
    canonicalAtomsWithOxt.insert("OXT");
    for (const auto& atom : sourceAtoms) {
        if (canonicalAtomsWithOxt.find(atom) == canonicalAtomsWithOxt.end()) {
            return true;
        }
    }
    return sourceAtoms != canonicalAtoms && sourceAtoms != canonicalAtomsWithOxt;
}

bool regionNeedsRawFallbackForTrace(const tcb::span<AtomCoordinate>& atoms) {
    std::vector<std::pair<size_t, size_t>> residues = splitResidueRanges(atoms);
    size_t residuesWithOxt = 0;
    for (size_t i = 0; i < residues.size(); i++) {
        bool hasOxt = false;
        for (size_t j = residues[i].first; j < residues[i].second; j++) {
            if (atoms[j].atom == "OXT") {
                hasOxt = true;
                break;
            }
        }
        if (!hasOxt) {
            continue;
        }
        residuesWithOxt++;
        if (i + 1 != residues.size()) {
            return true;
        }
    }
    if (residuesWithOxt > 1) {
        return true;
    }
    for (const auto& residue : residues) {
        if (residueNeedsRawFallbackForTrace(
                tcb::span<const AtomCoordinate>(atoms.data() + residue.first, residue.second - residue.first))) {
            return true;
        }
    }
    return false;
}

bool shouldStoreSmallMixedFragmentAsRawForTrace(
    const tcb::span<AtomCoordinate>& chainSpan,
    const std::vector<BackboneRegion>& regions
) {
    if (regions.size() <= 1) {
        return false;
    }
    if (splitResidueRanges(chainSpan).size() > 24) {
        return false;
    }
    for (const auto& region : regions) {
        tcb::span<AtomCoordinate> regionSpan(
            chainSpan.data() + region.start,
            chainSpan.data() + region.end
        );
        if (!region.encodable || regionNeedsRawFallbackForTrace(regionSpan)) {
            return true;
        }
    }
    return false;
}

bool decodeEntryToAtoms(const char* dataBuffer, size_t size, std::vector<AtomCoordinate>& atoms) {
    std::string title;
    return decodeStructureToAtoms(dataBuffer, size, false, title, atoms);
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

template <typename T>
bool readBinarySpan(const char* data, size_t size, size_t& offset, T* out, size_t count) {
    size_t bytes = sizeof(T) * count;
    if (offset + bytes > size) {
        return false;
    }
    memcpy(out, data + offset, bytes);
    offset += bytes;
    return true;
}

bool parseFczPayload(const std::string& payload, ParsedFczPayload& parsed) {
    parsed = ParsedFczPayload{};
    if (payload.size() < MAGICNUMBER_LENGTH + sizeof(CompressedFileHeader)) {
        return false;
    }
    if (memcmp(payload.data(), MAGICNUMBER, MAGICNUMBER_LENGTH) != 0) {
        return false;
    }
    size_t offset = MAGICNUMBER_LENGTH;
    if (!readBinaryValue(payload.data(), payload.size(), offset, parsed.header)) {
        return false;
    }

    parsed.anchorIndices.resize(parsed.header.nAnchor);
    if (!parsed.anchorIndices.empty() &&
        !readBinarySpan(payload.data(), payload.size(), offset,
                        parsed.anchorIndices.data(), parsed.anchorIndices.size())) {
        return false;
    }

    parsed.title.resize(parsed.header.lenTitle);
    if (!parsed.title.empty() &&
        !readBinarySpan(payload.data(), payload.size(), offset,
                        parsed.title.data(), parsed.title.size())) {
        return false;
    }

    parsed.prevCoords.resize(3);
    if (!readBinarySpan(payload.data(), payload.size(), offset,
                        parsed.prevCoords.data(), parsed.prevCoords.size())) {
        return false;
    }

    size_t innerAnchorCount = parsed.header.nAnchor > 2 ? (parsed.header.nAnchor - 2) * 3 : 0;
    parsed.innerAnchorCoords.resize(innerAnchorCount);
    if (!parsed.innerAnchorCoords.empty() &&
        !readBinarySpan(payload.data(), payload.size(), offset,
                        parsed.innerAnchorCoords.data(), parsed.innerAnchorCoords.size())) {
        return false;
    }

    parsed.lastCoords.resize(3);
    if (!readBinarySpan(payload.data(), payload.size(), offset,
                        parsed.lastCoords.data(), parsed.lastCoords.size())) {
        return false;
    }

    if (!readBinaryValue(payload.data(), payload.size(), offset, parsed.hasOXT) ||
        !readBinaryValue(payload.data(), payload.size(), offset, parsed.oxtCoords)) {
        return false;
    }

    parsed.compressedBackbone.resize(parsed.header.nResidue);
    if (!parsed.compressedBackbone.empty() &&
        !readBinarySpan(payload.data(), payload.size(), offset,
                        parsed.compressedBackbone.data(), parsed.compressedBackbone.size())) {
        return false;
    }

    parsed.sideChainAngles.resize(parsed.header.nSideChainTorsion);
    if (!parsed.sideChainAngles.empty() &&
        !readBinarySpan(payload.data(), payload.size(), offset,
                        parsed.sideChainAngles.data(), parsed.sideChainAngles.size())) {
        return false;
    }

    if (!readBinaryValue(payload.data(), payload.size(), offset, parsed.tempMin) ||
        !readBinaryValue(payload.data(), payload.size(), offset, parsed.tempContF)) {
        return false;
    }

    parsed.tempFactors.resize(parsed.header.nResidue);
    if (!parsed.tempFactors.empty() &&
        !readBinarySpan(payload.data(), payload.size(), offset,
                        parsed.tempFactors.data(), parsed.tempFactors.size())) {
        return false;
    }
    return offset == payload.size();
}

template <typename T>
std::string toStringValue(const T& value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

std::string formatFloatTrace(float value) {
    std::ostringstream oss;
    uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "float size mismatch");
    memcpy(&bits, &value, sizeof(bits));
    oss << std::setprecision(9) << value
        << " bits=0x" << std::hex << std::setw(8) << std::setfill('0') << bits
        << std::dec << std::setfill(' ');
    return oss.str();
}

bool firstHeaderFieldDiff(const ParsedFczPayload& left, const ParsedFczPayload& right, HeaderFieldDiff& diff) {
    auto setDiff = [&](const std::string& name, const auto& a, const auto& b) -> bool {
        if (a == b) {
            return false;
        }
        diff.name = name;
        diff.leftValue = toStringValue(a);
        diff.rightValue = toStringValue(b);
        return true;
    };

    if (setDiff("nResidue", left.header.nResidue, right.header.nResidue)) return true;
    if (setDiff("nAtom", left.header.nAtom, right.header.nAtom)) return true;
    if (setDiff("idxResidue", left.header.idxResidue, right.header.idxResidue)) return true;
    if (setDiff("idxAtom", left.header.idxAtom, right.header.idxAtom)) return true;
    if (setDiff("nAnchor", static_cast<int>(left.header.nAnchor), static_cast<int>(right.header.nAnchor))) return true;
    if (setDiff("chain", static_cast<int>(left.header.chain), static_cast<int>(right.header.chain))) return true;
    if (setDiff("version", static_cast<int>(left.header.version), static_cast<int>(right.header.version))) return true;
    if (setDiff("flags", static_cast<int>(left.header.flags), static_cast<int>(right.header.flags))) return true;
    if (setDiff("nSideChainTorsion", left.header.nSideChainTorsion, right.header.nSideChainTorsion)) return true;
    if (setDiff("firstResidue", static_cast<int>(left.header.firstResidue), static_cast<int>(right.header.firstResidue))) return true;
    if (setDiff("lastResidue", static_cast<int>(left.header.lastResidue), static_cast<int>(right.header.lastResidue))) return true;
    if (setDiff("chain2", static_cast<int>(left.header.chain2), static_cast<int>(right.header.chain2))) return true;
    if (setDiff("chain3", static_cast<int>(left.header.chain3), static_cast<int>(right.header.chain3))) return true;
    if (setDiff("lenTitle", left.header.lenTitle, right.header.lenTitle)) return true;
    for (size_t i = 0; i < 6; i++) {
        if (setDiff(std::string("mins[") + std::to_string(i) + "]",
                    left.header.mins[i], right.header.mins[i])) return true;
    }
    for (size_t i = 0; i < 6; i++) {
        if (setDiff(std::string("cont_fs[") + std::to_string(i) + "]",
                    left.header.cont_fs[i], right.header.cont_fs[i])) return true;
    }
    return false;
}

std::string sourcePathToDbKey(const std::string& path) {
    std::string base = baseName(path);
    std::pair<std::string, std::string> outputParts = getFileParts(base);
    std::string dbKey = baseName(outputParts.first);
    std::replace(dbKey.begin(), dbKey.end(), '.', '_');
    return dbKey;
}

std::vector<RawReason> classifyResidue(
    const std::vector<AtomCoordinate>& residueAtoms,
    bool isolatedCompleteResidue
) {
    std::vector<RawReason> reasons;
    std::set<std::string> atomNames;
    for (const auto& atom : residueAtoms) {
        atomNames.insert(atom.atom);
    }

    bool hasN = atomNames.count("N") != 0;
    bool hasCA = atomNames.count("CA") != 0;
    bool hasC = atomNames.count("C") != 0;
    if (atomNames.size() == 1 && hasCA) {
        reasons.push_back(RawReason::CaOnlyResidue);
    } else if (!(hasN && hasCA && hasC)) {
        reasons.push_back(RawReason::IncompleteBackboneResidue);
    }
    if (isolatedCompleteResidue) {
        reasons.push_back(RawReason::IsolatedCompleteResidue);
    }

    auto aaIt = Foldcomp::AAS.find(residueAtoms[0].residue);
    if (aaIt == Foldcomp::AAS.end()) {
        reasons.push_back(RawReason::UnknownResidue);
        return reasons;
    }
    std::set<std::string> canonicalAtoms(aaIt->second.atoms.begin(), aaIt->second.atoms.end());
    std::set<std::string> canonicalAtomsWithOxt = canonicalAtoms;
    canonicalAtomsWithOxt.insert("OXT");
    bool hasExtraNonCanonicalAtoms = false;
    for (const auto& atom : atomNames) {
        if (canonicalAtomsWithOxt.count(atom) == 0) {
            hasExtraNonCanonicalAtoms = true;
        }
    }
    if (hasExtraNonCanonicalAtoms) {
        reasons.push_back(RawReason::ExtraNonCanonicalAtoms);
    } else if (atomNames != canonicalAtoms && atomNames != canonicalAtomsWithOxt) {
        reasons.push_back(RawReason::PartialCanonicalAtomSet);
    }
    return reasons;
}

std::map<ResidueKey, std::vector<RawReason>> classifyRawResidues(const std::vector<AtomCoordinate>& atoms) {
    std::map<ResidueKey, std::vector<RawReason>> output;
    if (atoms.empty()) {
        return output;
    }
    std::map<std::tuple<int, std::string, int>, int> occurrences;

    std::vector<std::pair<size_t, size_t>> chainRegions = identifyChains(atoms);
    for (const auto& chainRegion : chainRegions) {
        std::vector<std::pair<size_t, size_t>> discontinuous = identifyDiscontinousResInd(
            atoms, chainRegion.first, chainRegion.second
        );
        if (discontinuous.empty()) {
            discontinuous.push_back(chainRegion);
        }
        for (const auto& fragment : discontinuous) {
            tcb::span<AtomCoordinate> fragmentSpan(
                const_cast<AtomCoordinate*>(&atoms[fragment.first]),
                const_cast<AtomCoordinate*>(atoms.data()) + fragment.second
            );
            std::vector<BackboneRegion> regions = identifyBackboneRegions(fragmentSpan);
            for (const auto& region : regions) {
                if (region.encodable) {
                    bool needsRaw = false;
                    std::vector<std::vector<AtomCoordinate>> residues = splitAtomByResidue(
                        tcb::span<AtomCoordinate>(fragmentSpan.data() + region.start,
                                                  fragmentSpan.data() + region.end)
                    );
                    for (const auto& residue : residues) {
                        auto reasons = classifyResidue(residue, false);
                        if (!reasons.empty()) {
                            needsRaw = true;
                            break;
                        }
                    }
                    if (!needsRaw) {
                        continue;
                    }
                }

                std::vector<std::vector<AtomCoordinate>> residues = splitAtomByResidue(
                    tcb::span<AtomCoordinate>(fragmentSpan.data() + region.start,
                                              fragmentSpan.data() + region.end)
                );
                for (const auto& residue : residues) {
                    bool atomicallyComplete = false;
                    {
                        std::set<std::string> names;
                        for (const auto& atom : residue) {
                            names.insert(atom.atom);
                        }
                        atomicallyComplete = names.count("N") && names.count("CA") && names.count("C");
                    }
                    bool isolatedComplete = region.encodable == false && atomicallyComplete;
                    auto reasons = classifyResidue(residue, isolatedComplete);
                    if (reasons.empty()) {
                        reasons.push_back(RawReason::PartialCanonicalAtomSet);
                    }
                    std::tuple<int, std::string, int> baseKey = std::make_tuple(
                        residue[0].model, residue[0].chain, residue[0].residue_index
                    );
                    ResidueKey key = std::make_tuple(
                        residue[0].model, residue[0].chain, residue[0].residue_index,
                        occurrences[baseKey]++
                    );
                    output[key] = reasons;
                }
            }
        }
    }
    return output;
}

int commandFallbackReport(const std::string& input) {
    FallbackStats stats;
    std::vector<std::string> files;
    if (fs::is_directory(input)) {
        for (const auto& entry : fs::directory_iterator(input)) {
            if (entry.is_regular_file()) {
                files.push_back(entry.path().string());
            }
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(input);
    }

    for (const auto& file : files) {
        stats.fileCount++;
        StructureReader reader;
        if (!reader.load(file)) {
            continue;
        }
        std::vector<AtomCoordinate> atoms;
        reader.readAllAtoms(atoms);
        removeAlternativePosition(atoms);
        if (atoms.empty()) {
            stats.rawFilesByReason[RawReason::NonProtein]++;
            continue;
        }
        stats.proteinFileCount++;
        auto rawResidues = classifyRawResidues(atoms);
        auto groupedResidues = groupResidues(atoms);
        stats.residueCount += groupedResidues.size();
        if (!rawResidues.empty()) {
            stats.rawFileCount++;
        }
        std::set<RawReason> fileReasons;
        for (const auto& kv : rawResidues) {
            for (RawReason reason : kv.second) {
                stats.rawResiduesByReason[reason]++;
                fileReasons.insert(reason);
            }
        }
        for (RawReason reason : fileReasons) {
            stats.rawFilesByReason[reason]++;
        }
    }

    std::cout << "files_total=" << stats.fileCount
              << " protein_files=" << stats.proteinFileCount
              << " raw_files=" << stats.rawFileCount
              << " residues_total=" << stats.residueCount
              << "\n";
    for (RawReason reason : {
             RawReason::NonProtein,
             RawReason::CaOnlyResidue,
             RawReason::IncompleteBackboneResidue,
             RawReason::IsolatedCompleteResidue,
             RawReason::UnknownResidue,
             RawReason::PartialCanonicalAtomSet,
             RawReason::ExtraNonCanonicalAtoms
         }) {
        size_t rawResidues = stats.rawResiduesByReason[reason];
        size_t rawFiles = stats.rawFilesByReason[reason];
        double residueRatio = stats.residueCount == 0 ? 0.0 : static_cast<double>(rawResidues) / static_cast<double>(stats.residueCount);
        double fileRatio = stats.proteinFileCount == 0 ? 0.0 : static_cast<double>(rawFiles) / static_cast<double>(stats.proteinFileCount);
        std::cout << rawReasonName(reason)
                  << " raw_files=" << rawFiles
                  << " file_ratio=" << fileRatio
                  << " raw_residues=" << rawResidues
                  << " residue_ratio=" << residueRatio
                  << "\n";
    }
    return 0;
}

Summary summarize(const std::vector<AtomCoordinate>& atoms) {
    std::map<ResidueKey, ResidueStats> residues = groupResidues(atoms);

    Summary summary;
    summary.atomCount = atoms.size();
    summary.residueCount = residues.size();
    for (const auto& kv : residues) {
        const auto& residue = kv.second;
        bool hasN = residue.atoms.count("N") != 0;
        bool hasCA = residue.atoms.count("CA") != 0;
        bool hasC = residue.atoms.count("C") != 0;
        if (hasN && hasCA && hasC) {
            summary.completeBackboneResidueCount++;
        }
        if (residue.atoms.size() == 1 && hasCA) {
            summary.caOnlyResidueCount++;
        }
    }
    return summary;
}

void printSummary(const std::string& label, const Summary& summary) {
    std::cout << label
              << " atoms=" << summary.atomCount
              << " residues=" << summary.residueCount
              << " complete_backbone=" << summary.completeBackboneResidueCount
              << " ca_only=" << summary.caOnlyResidueCount
              << "\n";
}

std::map<ResidueKey, ResidueStats> groupResidues(const std::vector<AtomCoordinate>& atoms) {
    std::map<ResidueKey, ResidueStats> residues;
    if (atoms.empty()) {
        return residues;
    }

    std::map<std::tuple<int, std::string, int, char, std::string>, int> occurrences;
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) ||
                            startsNewResidue(atoms[i], atoms[i - 1]);
        if (!endOfResidue) {
            continue;
        }

        const AtomCoordinate& firstAtom = atoms[residueStart];
        std::tuple<int, std::string, int, char, std::string> baseKey = std::make_tuple(
            firstAtom.model, firstAtom.chain, firstAtom.residue_index,
            firstAtom.insertion_code, firstAtom.residue
        );
        ResidueKey key = std::make_tuple(
            firstAtom.model, firstAtom.chain, firstAtom.residue_index,
            occurrences[baseKey]++
        );
        ResidueStats& stats = residues[key];
        stats.residue = firstAtom.residue;
        for (size_t j = residueStart; j < i; j++) {
            stats.atoms.insert(atoms[j].atom);
        }
        residueStart = i;
    }
    return residues;
}

std::string formatResidueKey(const ResidueKey& key, const ResidueStats& residue) {
    std::ostringstream oss;
    oss << "model=" << std::get<0>(key)
        << " chain=" << std::get<1>(key)
        << " resid=" << std::get<2>(key)
        << " occ=" << std::get<3>(key)
        << " resn=" << residue.residue;
    return oss.str();
}

std::string joinAtomNames(const std::set<std::string>& atoms) {
    std::ostringstream oss;
    bool first = true;
    for (const auto& atom : atoms) {
        if (!first) {
            oss << ",";
        }
        oss << atom;
        first = false;
    }
    return oss.str();
}

void printResidueStream(const std::vector<AtomCoordinate>& atoms) {
    if (atoms.empty()) {
        return;
    }
    std::map<std::tuple<int, std::string, int>, int> occurrences;
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) ||
                            atoms[i].model != atoms[i - 1].model ||
                            atoms[i].chain != atoms[i - 1].chain ||
                            atoms[i].residue_index != atoms[i - 1].residue_index;
        if (!endOfResidue) {
            continue;
        }
        const AtomCoordinate& atom = atoms[residueStart];
        std::tuple<int, std::string, int> baseKey = std::make_tuple(atom.model, atom.chain, atom.residue_index);
        int occ = occurrences[baseKey]++;
        std::cout << "residue model=" << atom.model
                  << " chain=" << atom.chain
                  << " resid=" << atom.residue_index
                  << " occ=" << occ
                  << " resn=" << atom.residue
                  << "\n";
        residueStart = i;
    }
}

int commandSummarize(const std::string& input) {
    StructureReader reader;
    if (!reader.load(input)) {
        std::cerr << "[Error] Failed to read " << input << "\n";
        return 1;
    }
    std::vector<AtomCoordinate> atoms;
    reader.readAllAtoms(atoms);
    removeAlternativePosition(atoms);
    printSummary(input, summarize(atoms));
    return 0;
}

int compareAtomVectors(const std::vector<AtomCoordinate>& sourceAtoms,
                       const std::vector<AtomCoordinate>& roundtripAtoms) {
    Summary sourceSummary = summarize(sourceAtoms);
    Summary roundtripSummary = summarize(roundtripAtoms);
    printSummary("source", sourceSummary);
    printSummary("roundtrip", roundtripSummary);

    auto sourceResidues = groupResidues(sourceAtoms);
    auto roundtripResidues = groupResidues(roundtripAtoms);

    size_t missingResidues = 0;
    for (const auto& kv : sourceResidues) {
        if (roundtripResidues.find(kv.first) == roundtripResidues.end()) {
            std::cout << "missing_residue " << formatResidueKey(kv.first, kv.second) << "\n";
            missingResidues++;
        }
    }
    size_t extraResidues = 0;
    for (const auto& kv : roundtripResidues) {
        if (sourceResidues.find(kv.first) == sourceResidues.end()) {
            std::cout << "extra_residue " << formatResidueKey(kv.first, kv.second) << "\n";
            extraResidues++;
        }
    }

    size_t changedResidues = 0;
    for (const auto& kv : sourceResidues) {
        auto it = roundtripResidues.find(kv.first);
        if (it == roundtripResidues.end()) {
            continue;
        }
        if (kv.second.atoms != it->second.atoms) {
            std::set<std::string> missingAtoms;
            std::set<std::string> extraAtoms;
            std::set_difference(kv.second.atoms.begin(), kv.second.atoms.end(),
                                it->second.atoms.begin(), it->second.atoms.end(),
                                std::inserter(missingAtoms, missingAtoms.begin()));
            std::set_difference(it->second.atoms.begin(), it->second.atoms.end(),
                                kv.second.atoms.begin(), kv.second.atoms.end(),
                                std::inserter(extraAtoms, extraAtoms.begin()));
            std::cout << "changed_residue " << formatResidueKey(kv.first, kv.second)
                      << " source_atoms=" << kv.second.atoms.size()
                      << " roundtrip_atoms=" << it->second.atoms.size()
                      << " roundtrip_resn=" << it->second.residue
                      << " missing_atom_names=" << joinAtomNames(missingAtoms)
                      << " extra_atom_names=" << joinAtomNames(extraAtoms)
                      << "\n";
            changedResidues++;
        }
    }

    std::cout << "missing_residue_count=" << missingResidues
              << " extra_residue_count=" << extraResidues
              << " changed_residue_count=" << changedResidues
              << "\n";
    return 0;
}

std::map<RmsdResidueKey, std::vector<AtomCoordinate>> groupResiduesForRmsd(
    const std::vector<AtomCoordinate>& atoms
) {
    std::map<RmsdResidueKey, std::vector<AtomCoordinate>> residues;
    if (atoms.empty()) {
        return residues;
    }
    std::map<std::tuple<int, std::string, int, char, std::string>, int> occurrences;
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) ||
                            startsNewResidue(atoms[i], atoms[i - 1]);
        if (!endOfResidue) {
            continue;
        }
        const AtomCoordinate& first = atoms[residueStart];
        std::tuple<int, std::string, int, char, std::string> baseKey(
            first.model, first.chain, first.residue_index, first.insertion_code, first.residue
        );
        RmsdResidueKey key(
            first.model, first.chain, first.residue_index, occurrences[baseKey]++
        );
        std::vector<AtomCoordinate>& residue = residues[key];
        residue.insert(residue.end(), atoms.begin() + residueStart, atoms.begin() + i);
        residueStart = i;
    }
    return residues;
}

bool calculateAlignedRmsd(
    const std::vector<AtomCoordinate>& atoms1,
    const std::vector<AtomCoordinate>& atoms2,
    float& backboneRmsd,
    float& allAtomRmsd
) {
    backboneRmsd = 0.0f;
    allAtomRmsd = 0.0f;
    if (atoms1.empty() || atoms2.empty() || atoms1.size() != atoms2.size()) {
        return false;
    }

    std::vector<AtomCoordinate> alignedAtoms1;
    std::vector<AtomCoordinate> alignedAtoms2;
    std::vector<AtomCoordinate> alignedBackbone1;
    std::vector<AtomCoordinate> alignedBackbone2;

    std::map<RmsdResidueKey, std::vector<AtomCoordinate>> residues1 = groupResiduesForRmsd(atoms1);
    std::map<RmsdResidueKey, std::vector<AtomCoordinate>> residues2 = groupResiduesForRmsd(atoms2);
    if (residues1.size() != residues2.size()) {
        return false;
    }

    auto atomSorter = [](const AtomCoordinate& a, const AtomCoordinate& b) {
        if (a.atom != b.atom) {
            return a.atom < b.atom;
        }
        if (a.residue != b.residue) {
            return a.residue < b.residue;
        }
        if (a.chain != b.chain) {
            return a.chain < b.chain;
        }
        return a.atom_index < b.atom_index;
    };

    for (const auto& kv : residues1) {
        auto it = residues2.find(kv.first);
        if (it == residues2.end() || kv.second.empty() || it->second.empty()) {
            return false;
        }
        const AtomCoordinate& left = kv.second.front();
        const AtomCoordinate& right = it->second.front();
        if (left.residue != right.residue || kv.second.size() != it->second.size()) {
            return false;
        }

        std::vector<AtomCoordinate> sorted1 = kv.second;
        std::vector<AtomCoordinate> sorted2 = it->second;
        std::sort(sorted1.begin(), sorted1.end(), atomSorter);
        std::sort(sorted2.begin(), sorted2.end(), atomSorter);
        for (size_t atomIdx = 0; atomIdx < sorted1.size(); atomIdx++) {
            if (sorted1[atomIdx].atom != sorted2[atomIdx].atom) {
                return false;
            }
            alignedAtoms1.push_back(sorted1[atomIdx]);
            alignedAtoms2.push_back(sorted2[atomIdx]);
        }

        auto appendBackboneAtom = [&](const std::string& atomName) {
            auto leftIt = std::find_if(
                sorted1.begin(), sorted1.end(),
                [&](const AtomCoordinate& atom) { return atom.atom == atomName; }
            );
            auto rightIt = std::find_if(
                sorted2.begin(), sorted2.end(),
                [&](const AtomCoordinate& atom) { return atom.atom == atomName; }
            );
            if (leftIt != sorted1.end() && rightIt != sorted2.end()) {
                alignedBackbone1.push_back(*leftIt);
                alignedBackbone2.push_back(*rightIt);
            }
        };
        appendBackboneAtom("N");
        appendBackboneAtom("CA");
        appendBackboneAtom("C");
    }

    if (alignedAtoms1.empty() || alignedBackbone1.empty()) {
        return false;
    }
    backboneRmsd = RMSD(alignedBackbone1, alignedBackbone2);
    allAtomRmsd = RMSD(alignedAtoms1, alignedAtoms2);
    return true;
}

int commandCompare(const std::string& source, const std::string& roundtrip) {
    StructureReader sourceReader;
    StructureReader roundtripReader;
    if (!sourceReader.load(source)) {
        std::cerr << "[Error] Failed to read source " << source << "\n";
        return 1;
    }
    if (!roundtripReader.load(roundtrip)) {
        std::cerr << "[Error] Failed to read roundtrip " << roundtrip << "\n";
        return 1;
    }

    std::vector<AtomCoordinate> sourceAtoms;
    std::vector<AtomCoordinate> roundtripAtoms;
    sourceReader.readAllAtoms(sourceAtoms);
    roundtripReader.readAllAtoms(roundtripAtoms);
    removeAlternativePosition(sourceAtoms);
    removeAlternativePosition(roundtripAtoms);
    return compareAtomVectors(sourceAtoms, roundtripAtoms);
}

int commandDbEntryCompare(const std::string& sourcePath, const std::string& dbPath) {
    StructureReader sourceReader;
    if (!sourceReader.load(sourcePath)) {
        std::cerr << "[Error] Failed to read source " << sourcePath << "\n";
        return 1;
    }

    std::vector<AtomCoordinate> sourceAtoms;
    sourceReader.readAllAtoms(sourceAtoms);
    removeAlternativePosition(sourceAtoms);
    if (sourceAtoms.empty()) {
        std::cerr << "[Error] No protein atoms found in source " << sourcePath << "\n";
        return 1;
    }

    void* reader = make_reader(dbPath.c_str(), (dbPath + ".index").c_str(),
                               DB_READER_USE_DATA | DB_READER_USE_LOOKUP);
    if (reader == NULL) {
        std::cerr << "[Error] Failed to open database " << dbPath << "\n";
        return 1;
    }

    std::string dbKeyStr = sourcePathToDbKey(sourcePath);
    uint32_t key = reader_lookup_entry(reader, dbKeyStr.c_str());
    if (key == UINT32_MAX) {
        std::cerr << "[Error] Missing db entry for key " << dbKeyStr << "\n";
        free_reader(reader);
        return 1;
    }
    int64_t id = reader_get_id(reader, key);
    if (id == -1) {
        std::cerr << "[Error] Missing db id for key " << dbKeyStr << "\n";
        free_reader(reader);
        return 1;
    }

    std::vector<AtomCoordinate> dbAtoms;
    const char* data = reader_get_data(reader, id);
    int64_t length = reader_get_length(reader, id);
    if (data == NULL || length <= 0 || !decodeEntryToAtoms(data, static_cast<size_t>(length), dbAtoms)) {
        std::cerr << "[Error] Failed to decode db entry for key " << dbKeyStr << "\n";
        free_reader(reader);
        return 1;
    }
    free_reader(reader);

    removeAlternativePosition(dbAtoms);
    return compareAtomVectors(sourceAtoms, dbAtoms);
}

std::string bucketizeRmsd(float value, const std::vector<float>& bins) {
    float lower = bins.front();
    for (size_t i = 1; i < bins.size(); i++) {
        float upper = bins[i];
        if (value <= upper) {
            std::ostringstream oss;
            if (std::isinf(upper)) {
                oss << ">" << std::fixed << std::setprecision(2) << lower;
            } else {
                oss << std::fixed << std::setprecision(2) << lower
                    << "-" << std::fixed << std::setprecision(2) << upper;
            }
            return oss.str();
        }
        lower = upper;
    }
    return "unreachable";
}

void printHistogram(const std::string& label, const std::map<std::string, size_t>& histogram) {
    for (const auto& kv : histogram) {
        std::cout << label << " bucket=" << kv.first << " count=" << kv.second << "\n";
    }
}

int commandDbRmsd(const std::string& sourceDir, const std::string& dbPath, bool recursive) {
    std::vector<std::string> files = getFilesInDirectory(sourceDir, recursive);
    std::sort(files.begin(), files.end());

    void* reader = make_reader(dbPath.c_str(), (dbPath + ".index").c_str(),
                               DB_READER_USE_DATA | DB_READER_USE_LOOKUP);
    if (reader == NULL) {
        std::cerr << "[Error] Failed to open database " << dbPath << "\n";
        return 1;
    }

    const std::vector<float> backboneBins = {0.0f, 0.01f, 0.02f, 0.05f, 0.1f, 0.2f, 0.5f, 1.0f, std::numeric_limits<float>::infinity()};
    const std::vector<float> allAtomBins = {0.0f, 0.02f, 0.05f, 0.1f, 0.2f, 0.5f, 1.0f, 2.0f, std::numeric_limits<float>::infinity()};

    size_t totalFiles = 0;
    size_t proteinFiles = 0;
    size_t matchedFiles = 0;
    size_t missingEntries = 0;
    size_t decodeFailures = 0;
    size_t rmsdFailures = 0;
    std::vector<float> backboneValues;
    std::vector<float> allAtomValues;
    std::map<std::string, size_t> backboneHistogram;
    std::map<std::string, size_t> allAtomHistogram;

    for (const auto& file : files) {
        totalFiles++;
        StructureReader sourceReader;
        if (!sourceReader.load(file)) {
            continue;
        }
        std::vector<AtomCoordinate> sourceAtoms;
        sourceReader.readAllAtoms(sourceAtoms);
        removeAlternativePosition(sourceAtoms);
        if (sourceAtoms.empty()) {
            continue;
        }
        proteinFiles++;

        std::string dbKeyStr = sourcePathToDbKey(file);
        uint32_t key = reader_lookup_entry(reader, dbKeyStr.c_str());
        if (key == UINT32_MAX) {
            std::cout << "missing_db_entry " << baseName(file) << " key=" << dbKeyStr << "\n";
            missingEntries++;
            continue;
        }
        int64_t id = reader_get_id(reader, key);
        if (id == -1) {
            std::cout << "missing_db_id " << baseName(file) << " key=" << dbKeyStr << "\n";
            missingEntries++;
            continue;
        }

        std::vector<AtomCoordinate> dbAtoms;
        const char* data = reader_get_data(reader, id);
        int64_t length = reader_get_length(reader, id);
        if (data == NULL || length <= 0 || !decodeEntryToAtoms(data, static_cast<size_t>(length), dbAtoms)) {
            std::cout << "decode_failure " << baseName(file) << " key=" << dbKeyStr << "\n";
            decodeFailures++;
            continue;
        }
        removeAlternativePosition(dbAtoms);

        float backboneRmsd = 0.0f;
        float allAtomRmsd = 0.0f;
        if (!calculateAlignedRmsd(sourceAtoms, dbAtoms, backboneRmsd, allAtomRmsd)) {
            std::cout << "rmsd_failure " << baseName(file) << " key=" << dbKeyStr << "\n";
            rmsdFailures++;
            continue;
        }

        matchedFiles++;
        backboneValues.push_back(backboneRmsd);
        allAtomValues.push_back(allAtomRmsd);
        backboneHistogram[bucketizeRmsd(backboneRmsd, backboneBins)]++;
        allAtomHistogram[bucketizeRmsd(allAtomRmsd, allAtomBins)]++;
    }
    free_reader(reader);

    auto printStats = [](const std::string& prefix, const std::vector<float>& values) {
        if (values.empty()) {
            return;
        }
        std::vector<float> sorted = values;
        std::sort(sorted.begin(), sorted.end());
        float mean = std::accumulate(sorted.begin(), sorted.end(), 0.0f) / static_cast<float>(sorted.size());
        float median = sorted[sorted.size() / 2];
        if ((sorted.size() % 2) == 0) {
            median = 0.5f * (sorted[sorted.size() / 2 - 1] + sorted[sorted.size() / 2]);
        }
        std::cout << prefix
                  << " mean=" << mean
                  << " median=" << median
                  << " max=" << sorted.back()
                  << "\n";
    };

    std::cout << "files_total=" << totalFiles
              << " protein_files=" << proteinFiles
              << " matched_files=" << matchedFiles
              << " missing_entries=" << missingEntries
              << " decode_failures=" << decodeFailures
              << " rmsd_failures=" << rmsdFailures
              << "\n";
    printStats("backbone_rmsd", backboneValues);
    printStats("all_atom_rmsd", allAtomValues);
    printHistogram("backbone_hist", backboneHistogram);
    printHistogram("all_atom_hist", allAtomHistogram);
    return 0;
}

int commandDumpResidues(const std::string& sourcePath) {
    StructureReader reader;
    if (!reader.load(sourcePath)) {
        std::cerr << "[Error] Failed to read source " << sourcePath << "\n";
        return 1;
    }
    std::vector<AtomCoordinate> atoms;
    reader.readAllAtoms(atoms);
    removeAlternativePosition(atoms);
    printResidueStream(atoms);
    return 0;
}

int commandDbEntryDumpResidues(const std::string& sourcePath, const std::string& dbPath) {
    void* reader = make_reader(dbPath.c_str(), (dbPath + ".index").c_str(),
                               DB_READER_USE_DATA | DB_READER_USE_LOOKUP);
    if (reader == NULL) {
        std::cerr << "[Error] Failed to open database " << dbPath << "\n";
        return 1;
    }
    std::string dbKeyStr = sourcePathToDbKey(sourcePath);
    uint32_t key = reader_lookup_entry(reader, dbKeyStr.c_str());
    int64_t id = (key == UINT32_MAX) ? -1 : reader_get_id(reader, key);
    std::vector<AtomCoordinate> atoms;
    const char* data = (id == -1) ? NULL : reader_get_data(reader, id);
    int64_t length = (id == -1) ? -1 : reader_get_length(reader, id);
    if (id == -1 || data == NULL || length <= 0 ||
        !decodeEntryToAtoms(data, static_cast<size_t>(length), atoms)) {
        std::cerr << "[Error] Failed to decode db entry for key " << dbKeyStr << "\n";
        free_reader(reader);
        return 1;
    }
    free_reader(reader);
    removeAlternativePosition(atoms);
    printResidueStream(atoms);
    return 0;
}

int commandFragmentReport(const std::string& sourcePath) {
    StructureReader reader;
    if (!reader.load(sourcePath)) {
        std::cerr << "[Error] Failed to read source " << sourcePath << "\n";
        return 1;
    }
    std::vector<AtomCoordinate> atoms;
    reader.readAllAtoms(atoms);
    removeAlternativePosition(atoms);
    if (atoms.empty()) {
        std::cerr << "[Error] No protein atoms found in " << sourcePath << "\n";
        return 1;
    }

    std::vector<std::pair<size_t, size_t>> chainRegions = identifyChains(atoms);
    for (size_t chainIdx = 0; chainIdx < chainRegions.size(); chainIdx++) {
        const auto& chainRegion = chainRegions[chainIdx];
        std::cout << "chain_region idx=" << chainIdx
                  << " start_atom=" << chainRegion.first
                  << " end_atom=" << chainRegion.second
                  << " start_resid=" << atoms[chainRegion.first].residue_index
                  << " end_resid=" << atoms[chainRegion.second - 1].residue_index
                  << " chain=" << atoms[chainRegion.first].chain
                  << "\n";
        std::vector<std::pair<size_t, size_t>> fragments = identifyDiscontinousResInd(
            atoms, chainRegion.first, chainRegion.second
        );
        if (fragments.empty()) {
            fragments.push_back(chainRegion);
        }
        for (size_t fragIdx = 0; fragIdx < fragments.size(); fragIdx++) {
            const auto& fragment = fragments[fragIdx];
            tcb::span<AtomCoordinate> chainSpan(
                &atoms[fragment.first], atoms.data() + fragment.second
            );
            std::vector<BackboneRegion> regions = identifyBackboneRegions(chainSpan);
            std::cout << "fragment idx=" << fragIdx
                      << " start_atom=" << fragment.first
                      << " end_atom=" << fragment.second
                      << " start_resid=" << atoms[fragment.first].residue_index
                      << " end_resid=" << atoms[fragment.second - 1].residue_index
                      << " regions=" << regions.size()
                      << "\n";
            for (size_t regionIdx = 0; regionIdx < regions.size(); regionIdx++) {
                const auto& region = regions[regionIdx];
                size_t startAtom = fragment.first + region.start;
                size_t endAtom = fragment.first + region.end - 1;
                std::cout << "region idx=" << regionIdx
                          << " encodable=" << (region.encodable ? 1 : 0)
                          << " start_atom=" << startAtom
                          << " end_atom=" << endAtom
                          << " start_resid=" << atoms[startAtom].residue_index
                          << " end_resid=" << atoms[endAtom].residue_index
                          << "\n";
            }
        }
    }
    return 0;
}

bool findTraceFragment(const std::vector<AtomCoordinate>& atoms, int selectedFragmentIndex, TraceFragmentSelection& selection) {
    int globalFragmentIndex = 0;
    std::vector<std::pair<size_t, size_t>> chainRegions = identifyChains(atoms);
    for (size_t chainIdx = 0; chainIdx < chainRegions.size(); chainIdx++) {
        const auto& chainRegion = chainRegions[chainIdx];
        std::vector<std::pair<size_t, size_t>> fragments = identifyDiscontinousResInd(
            atoms, chainRegion.first, chainRegion.second
        );
        if (fragments.empty()) {
            fragments.push_back(chainRegion);
        }
        for (size_t fragIdx = 0; fragIdx < fragments.size(); fragIdx++) {
            const auto& fragment = fragments[fragIdx];
            tcb::span<AtomCoordinate> chainSpan(
                const_cast<AtomCoordinate*>(&atoms[fragment.first]), fragment.second - fragment.first
            );
            std::vector<BackboneRegion> regions = identifyBackboneRegions(chainSpan);
            bool forceRawWholeChain = shouldStoreSmallMixedFragmentAsRawForTrace(chainSpan, regions);
            for (size_t regionIdx = 0; regionIdx < regions.size(); regionIdx++, globalFragmentIndex++) {
                if (globalFragmentIndex != selectedFragmentIndex) {
                    continue;
                }
                const auto& region = regions[regionIdx];
                tcb::span<AtomCoordinate> regionSpan(
                    chainSpan.data() + region.start,
                    chainSpan.data() + region.end
                );
                selection.globalFragmentIndex = globalFragmentIndex;
                selection.chainRegionIndex = chainIdx;
                selection.discontinuityIndex = fragIdx;
                selection.regionIndex = regionIdx;
                selection.useFoldcompEncoding = !forceRawWholeChain &&
                                               region.encodable &&
                                               !regionNeedsRawFallbackForTrace(regionSpan);
                selection.regionAtoms.assign(regionSpan.begin(), regionSpan.end());
                return true;
            }
        }
    }
    return false;
}

int commandTraceEncode(const std::string& sourcePath, int selectedFragmentIndex) {
    StructureReader reader;
    if (!reader.load(sourcePath)) {
        std::cerr << "[Error] Failed to read source " << sourcePath << "\n";
        return 1;
    }
    std::vector<AtomCoordinate> atoms;
    reader.readAllAtoms(atoms);
    if (atoms.empty()) {
        std::cerr << "[Error] No protein atoms found in " << sourcePath << "\n";
        return 1;
    }

    bool found = false;
    std::vector<std::pair<size_t, size_t>> chainRegions = identifyChains(atoms);
    int globalFragmentIndex = 0;
    for (size_t chainIdx = 0; chainIdx < chainRegions.size(); chainIdx++) {
        const auto& chainRegion = chainRegions[chainIdx];
        std::vector<std::pair<size_t, size_t>> fragments = identifyDiscontinousResInd(
            atoms, chainRegion.first, chainRegion.second
        );
        if (fragments.empty()) {
            fragments.push_back(chainRegion);
        }
        for (size_t fragIdx = 0; fragIdx < fragments.size(); fragIdx++) {
            const auto& fragment = fragments[fragIdx];
            tcb::span<AtomCoordinate> chainSpan(
                &atoms[fragment.first], atoms.data() + fragment.second
            );
            std::vector<BackboneRegion> regions = identifyBackboneRegions(chainSpan);
            bool forceRawWholeChain = shouldStoreSmallMixedFragmentAsRawForTrace(chainSpan, regions);
            for (size_t regionIdx = 0; regionIdx < regions.size(); regionIdx++, globalFragmentIndex++) {
                const auto& region = regions[regionIdx];
                tcb::span<AtomCoordinate> regionSpan(
                    chainSpan.data() + region.start,
                    chainSpan.data() + region.end
                );
                bool useFoldcompEncoding = !forceRawWholeChain &&
                                           region.encodable &&
                                           !regionNeedsRawFallbackForTrace(regionSpan);
                if (selectedFragmentIndex >= 0 && globalFragmentIndex != selectedFragmentIndex) {
                    continue;
                }
                found = true;
                size_t startAtom = fragment.first + region.start;
                size_t endAtom = fragment.first + region.end - 1;
                std::cout << "fragment idx=" << globalFragmentIndex
                          << " chain_region_idx=" << chainIdx
                          << " discontinuity_idx=" << fragIdx
                          << " region_idx=" << regionIdx
                          << " kind=" << (useFoldcompEncoding ? "fcz" : "raw")
                          << " chain=" << atoms[startAtom].chain
                          << " model=" << atoms[startAtom].model
                          << " start_atom=" << startAtom
                          << " end_atom=" << endAtom
                          << " start_resid=" << atoms[startAtom].residue_index
                          << " end_resid=" << atoms[endAtom].residue_index
                          << " atom_count=" << regionSpan.size()
                          << "\n";
                if (!useFoldcompEncoding) {
                    continue;
                }

                Foldcomp compRes;
                compRes.compress(regionSpan);
                std::cout << "header idxResidue=" << compRes.header.idxResidue
                          << " nResidue=" << compRes.header.nResidue
                          << " idxAtom=" << compRes.header.idxAtom
                          << " nAtom=" << compRes.header.nAtom
                          << " nAnchor=" << static_cast<int>(compRes.header.nAnchor)
                          << " nSideChainTorsion=" << compRes.header.nSideChainTorsion
                          << "\n";
                std::cout << "mins";
                for (int i = 0; i < 6; i++) {
                    std::cout << " [" << i << "]=" << formatFloatTrace(compRes.header.mins[i]);
                }
                std::cout << "\n";
                std::cout << "cont_fs";
                for (int i = 0; i < 6; i++) {
                    std::cout << " [" << i << "]=" << formatFloatTrace(compRes.header.cont_fs[i]);
                }
                std::cout << "\n";
                std::cout << "anchor_indices";
                for (int anchorIndex : compRes.anchorIndices) {
                    std::cout << " " << anchorIndex;
                }
                std::cout << "\n";

                for (size_t i = 0; i < compRes.compressedBackBone.size(); i++) {
                    const BackboneChain& bb = compRes.compressedBackBone[i];
                    std::string residue3 = i < compRes.residueThreeLetter.size() ? compRes.residueThreeLetter[i] : "";
                    std::cout << "bb idx=" << i
                              << " residue=" << residue3
                              << " residue_disc=" << bb.residue;
                    if (i < compRes.phi.size()) {
                        std::cout << " phi=" << formatFloatTrace(compRes.phi[i])
                                  << " phi_disc=" << compRes.phiDiscretized[i]
                                  << " psi=" << formatFloatTrace(compRes.psi[i])
                                  << " psi_disc=" << compRes.psiDiscretized[i]
                                  << " omega=" << formatFloatTrace(compRes.omega[i])
                                  << " omega_disc=" << compRes.omegaDiscretized[i]
                                  << " n_ca_c=" << formatFloatTrace(compRes.n_ca_c_angle[i])
                                  << " n_ca_c_disc=" << compRes.n_ca_c_angleDiscretized[i]
                                  << " ca_c_n=" << formatFloatTrace(compRes.ca_c_n_angle[i])
                                  << " ca_c_n_disc=" << compRes.ca_c_n_angleDiscretized[i]
                                  << " c_n_ca=" << formatFloatTrace(compRes.c_n_ca_angle[i])
                                  << " c_n_ca_disc=" << compRes.c_n_ca_angleDiscretized[i];
                    } else {
                        std::cout << " terminal=1";
                    }
                    std::cout << "\n";
                }

                size_t sidechainFlatIndex = 0;
                for (size_t residueIndex = 0; residueIndex < compRes.sideChainAnglesPerResidue.size(); residueIndex++) {
                    for (size_t torsionIndex = 0; torsionIndex < compRes.sideChainAnglesPerResidue[residueIndex].size(); torsionIndex++, sidechainFlatIndex++) {
                        std::cout << "sc residue_idx=" << residueIndex
                                  << " torsion_idx=" << torsionIndex
                                  << " value=" << formatFloatTrace(compRes.sideChainAnglesPerResidue[residueIndex][torsionIndex]);
                        if (sidechainFlatIndex < compRes.sideChainAnglesDiscretized.size()) {
                            std::cout << " disc=" << compRes.sideChainAnglesDiscretized[sidechainFlatIndex];
                        }
                        std::cout << "\n";
                    }
                }

                for (size_t i = 0; i < compRes.tempFactors.size(); i++) {
                    std::cout << "tf idx=" << i
                              << " value=" << formatFloatTrace(compRes.tempFactors[i]);
                    if (i < compRes.tempFactorsDiscretized.size()) {
                        std::cout << " disc=" << compRes.tempFactorsDiscretized[i];
                    }
                    std::cout << "\n";
                }
            }
        }
    }
    if (!found) {
        std::cerr << "[Error] Fragment index " << selectedFragmentIndex << " not found\n";
        return 1;
    }
    return 0;
}

int commandTraceTorsionOp(const std::string& sourcePath, int selectedFragmentIndex, int torsionIndex) {
    StructureReader reader;
    if (!reader.load(sourcePath)) {
        std::cerr << "[Error] Failed to read source " << sourcePath << "\n";
        return 1;
    }
    std::vector<AtomCoordinate> atoms;
    reader.readAllAtoms(atoms);
    if (atoms.empty()) {
        std::cerr << "[Error] No protein atoms found in " << sourcePath << "\n";
        return 1;
    }

    TraceFragmentSelection selection;
    if (!findTraceFragment(atoms, selectedFragmentIndex, selection)) {
        std::cerr << "[Error] Fragment index " << selectedFragmentIndex << " not found\n";
        return 1;
    }
    if (!selection.useFoldcompEncoding) {
        std::cerr << "[Error] Fragment index " << selectedFragmentIndex << " is not FCZ-encodable\n";
        return 1;
    }

    Foldcomp compRes;
    compRes.compress(selection.regionAtoms);
    std::vector<float> torsionVector = getTorsionFromXYZ(compRes.backbone, 1);
    if (torsionIndex < 0 || static_cast<size_t>(torsionIndex + 3) >= compRes.backbone.size()) {
        std::cerr << "[Error] torsion_index out of range for fragment\n";
        return 1;
    }

    const float3d atm1 = compRes.backbone[torsionIndex + 0].coordinate;
    const float3d atm2 = compRes.backbone[torsionIndex + 1].coordinate;
    const float3d atm3 = compRes.backbone[torsionIndex + 2].coordinate;
    const float3d atm4 = compRes.backbone[torsionIndex + 3].coordinate;
    const float3d d1{atm2.x - atm1.x, atm2.y - atm1.y, atm2.z - atm1.z};
    const float3d d2{atm3.x - atm2.x, atm3.y - atm2.y, atm3.z - atm2.z};
    const float3d d3{atm4.x - atm3.x, atm4.y - atm3.y, atm4.z - atm3.z};
    const float3d u1 = crossProduct(d1, d2);
    const float3d u2 = crossProduct(d2, d3);
    const float d2Norm = norm(d2);
    const float x = dotProduct(u1, u2);
    const float3d planeBetaVec = crossProduct(u2, d2);
    const float y = d2Norm > FLOAT3D_EPSILON ? dotProduct(u1, planeBetaVec) / d2Norm : 0.0f;
    const TorsionTrace inlineTrace = computeTorsionTraceFromCoords(atm1, atm2, atm3, atm4);
    const TorsionTrace runtimeTrace = computeTorsionTraceRuntime(atm1, atm2, atm3, atm4);
    const float torsionAngle = inlineTrace.angle;
    const int angleKind = torsionIndex % 3;
    const int residueIndex = torsionIndex / 3 + (angleKind == 2 ? 1 : 0);
    const int storedIndex = torsionIndex / 3;
    const char* angleName = angleKind == 0 ? "psi" : (angleKind == 1 ? "omega" : "phi");

    auto printVec = [](const char* name, const float3d& value) {
        std::cout << name
                  << " x=" << formatFloatTrace(value.x)
                  << " y=" << formatFloatTrace(value.y)
                  << " z=" << formatFloatTrace(value.z)
                  << "\n";
    };

    std::cout << "fragment idx=" << selection.globalFragmentIndex
              << " torsion_index=" << torsionIndex
              << " angle=" << angleName
              << " residue_idx=" << residueIndex
              << "\n";
    printVec("atm1", atm1);
    printVec("atm2", atm2);
    printVec("atm3", atm3);
    printVec("atm4", atm4);
    printVec("d1", d1);
    printVec("d2", d2);
    printVec("d3", d3);
    printVec("u1", u1);
    printVec("u2", u2);
    printVec("plane_beta_vec", planeBetaVec);
    std::cout << "d2_norm=" << formatFloatTrace(d2Norm) << "\n";
    std::cout << "x=" << formatFloatTrace(x) << "\n";
    std::cout << "y=" << formatFloatTrace(y) << "\n";
    std::cout << "inline_trace d2_norm=" << formatFloatTrace(inlineTrace.d2Norm)
              << " x=" << formatFloatTrace(inlineTrace.x)
              << " y=" << formatFloatTrace(inlineTrace.y)
              << " angle=" << formatFloatTrace(inlineTrace.angle)
              << "\n";
    std::cout << "runtime_trace d2_norm=" << formatFloatTrace(runtimeTrace.d2Norm)
              << " x=" << formatFloatTrace(runtimeTrace.x)
              << " y=" << formatFloatTrace(runtimeTrace.y)
              << " angle=" << formatFloatTrace(runtimeTrace.angle)
              << "\n";
    if (angleKind == 0 && static_cast<size_t>(storedIndex) < compRes.psi.size()) {
        std::cout << "stored_angle=" << formatFloatTrace(compRes.psi[storedIndex]) << "\n";
    } else if (angleKind == 1 && static_cast<size_t>(storedIndex) < compRes.omega.size()) {
        std::cout << "stored_angle=" << formatFloatTrace(compRes.omega[storedIndex]) << "\n";
    } else if (angleKind == 2 && static_cast<size_t>(storedIndex) < compRes.phi.size()) {
        std::cout << "stored_angle=" << formatFloatTrace(compRes.phi[storedIndex]) << "\n";
    }
    if (static_cast<size_t>(torsionIndex) < torsionVector.size()) {
        std::cout << "function_angle=" << formatFloatTrace(torsionVector[torsionIndex]) << "\n";
    }
    std::cout << "torsion_angle=" << formatFloatTrace(torsionAngle) << "\n";
    return 0;
}

int commandEncodedDiff(const std::string& leftPath, const std::string& rightPath) {
    std::string leftData;
    std::string rightData;
    {
        std::ifstream input(leftPath, std::ios::binary);
        if (!input) {
            std::cerr << "[Error] Failed to read " << leftPath << "\n";
            return 1;
        }
        leftData.assign((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    }
    {
        std::ifstream input(rightPath, std::ios::binary);
        if (!input) {
            std::cerr << "[Error] Failed to read " << rightPath << "\n";
            return 1;
        }
        rightData.assign((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    }
    if (leftData.empty()) {
        std::cerr << "[Error] Failed to read " << leftPath << "\n";
        return 1;
    }
    if (rightData.empty()) {
        std::cerr << "[Error] Failed to read " << rightPath << "\n";
        return 1;
    }

    std::string leftTitle;
    std::string rightTitle;
    std::vector<ContainerFragment> leftFragments;
    std::vector<ContainerFragment> rightFragments;
    bool leftContainer = readContainer(leftData.data(), leftData.size(), leftTitle, leftFragments);
    bool rightContainer = readContainer(rightData.data(), rightData.size(), rightTitle, rightFragments);

    if (!leftContainer) {
        ContainerFragment fragment;
        fragment.kind = CONTAINER_FRAGMENT_KIND_FCZ;
        fragment.payload = leftData;
        leftFragments.push_back(std::move(fragment));
    }
    if (!rightContainer) {
        ContainerFragment fragment;
        fragment.kind = CONTAINER_FRAGMENT_KIND_FCZ;
        fragment.payload = rightData;
        rightFragments.push_back(std::move(fragment));
    }

    std::cout << "left_container=" << (leftContainer ? 1 : 0)
              << " right_container=" << (rightContainer ? 1 : 0)
              << " left_fragments=" << leftFragments.size()
              << " right_fragments=" << rightFragments.size()
              << "\n";

    size_t fragmentCount = std::min(leftFragments.size(), rightFragments.size());
    for (size_t i = 0; i < fragmentCount; i++) {
        const auto& leftFragment = leftFragments[i];
        const auto& rightFragment = rightFragments[i];
        bool samePayload = leftFragment.payload == rightFragment.payload;
        size_t firstByteDiff = 0;
        bool hasFirstByteDiff = false;
        size_t payloadCompareCount = std::min(leftFragment.payload.size(), rightFragment.payload.size());
        for (size_t j = 0; j < payloadCompareCount; j++) {
            if (leftFragment.payload[j] != rightFragment.payload[j]) {
                firstByteDiff = j;
                hasFirstByteDiff = true;
                break;
            }
        }
        std::cout << "fragment idx=" << i
                  << " left_kind=" << static_cast<int>(leftFragment.kind)
                  << " right_kind=" << static_cast<int>(rightFragment.kind)
                  << " left_model=" << leftFragment.model
                  << " right_model=" << rightFragment.model
                  << " left_chain=" << leftFragment.chain
                  << " right_chain=" << rightFragment.chain
                  << " same_payload=" << (samePayload ? 1 : 0)
                  << " left_size=" << leftFragment.payload.size()
                  << " right_size=" << rightFragment.payload.size();
        if (hasFirstByteDiff) {
            std::cout << " first_payload_byte_diff=" << firstByteDiff;
        }
        std::cout
                  << "\n";
        if (samePayload || leftFragment.kind != CONTAINER_FRAGMENT_KIND_FCZ ||
            rightFragment.kind != CONTAINER_FRAGMENT_KIND_FCZ) {
            continue;
        }

        ParsedFczPayload leftParsed;
        ParsedFczPayload rightParsed;
        if (!parseFczPayload(leftFragment.payload, leftParsed) ||
            !parseFczPayload(rightFragment.payload, rightParsed)) {
            std::cout << "fragment_parse_failed idx=" << i << "\n";
            continue;
        }

        std::cout << "fcz idx=" << i
                  << " left_nResidue=" << leftParsed.header.nResidue
                  << " right_nResidue=" << rightParsed.header.nResidue
                  << " left_nAtom=" << leftParsed.header.nAtom
                  << " right_nAtom=" << rightParsed.header.nAtom
                  << " left_nAnchor=" << static_cast<int>(leftParsed.header.nAnchor)
                  << " right_nAnchor=" << static_cast<int>(rightParsed.header.nAnchor)
                  << " left_nSideChain=" << leftParsed.header.nSideChainTorsion
                  << " right_nSideChain=" << rightParsed.header.nSideChainTorsion
                  << "\n";

        size_t headerByteDiffs = 0;
        const unsigned char* leftHeaderBytes = reinterpret_cast<const unsigned char*>(&leftParsed.header);
        const unsigned char* rightHeaderBytes = reinterpret_cast<const unsigned char*>(&rightParsed.header);
        for (size_t j = 0; j < sizeof(CompressedFileHeader); j++) {
            if (leftHeaderBytes[j] != rightHeaderBytes[j]) {
                headerByteDiffs++;
            }
        }
        if (headerByteDiffs > 0) {
            std::cout << "header_byte_diffs=" << headerByteDiffs << "\n";
            HeaderFieldDiff fieldDiff;
            if (firstHeaderFieldDiff(leftParsed, rightParsed, fieldDiff)) {
                std::cout << "first_header_field_diff name=" << fieldDiff.name
                          << " left=" << fieldDiff.leftValue
                          << " right=" << fieldDiff.rightValue
                          << "\n";
            }
        }

        size_t backboneDiffResidues = 0;
        size_t residueDiff = 0;
        size_t phiDiff = 0;
        size_t psiDiff = 0;
        size_t omegaDiff = 0;
        size_t nCaCDiff = 0;
        size_t caCNDiff = 0;
        size_t cNCaDiff = 0;
        std::vector<std::string> examples;
        size_t backboneCount = std::min(leftParsed.compressedBackbone.size(), rightParsed.compressedBackbone.size());
        for (size_t j = 0; j < backboneCount; j++) {
            const BackboneChain& leftBb = leftParsed.compressedBackbone[j];
            const BackboneChain& rightBb = rightParsed.compressedBackbone[j];
            bool differs = false;
            std::string firstFieldName;
            std::string firstFieldLeft;
            std::string firstFieldRight;
            if (leftBb.residue != rightBb.residue) {
                residueDiff++;
                differs = true;
                if (firstFieldName.empty()) {
                    firstFieldName = "residue";
                    firstFieldLeft = toStringValue(static_cast<int>(leftBb.residue));
                    firstFieldRight = toStringValue(static_cast<int>(rightBb.residue));
                }
            }
            if (leftBb.phi != rightBb.phi) {
                phiDiff++;
                differs = true;
                if (firstFieldName.empty()) {
                    firstFieldName = "phi";
                    firstFieldLeft = toStringValue(leftBb.phi);
                    firstFieldRight = toStringValue(rightBb.phi);
                }
            }
            if (leftBb.psi != rightBb.psi) {
                psiDiff++;
                differs = true;
                if (firstFieldName.empty()) {
                    firstFieldName = "psi";
                    firstFieldLeft = toStringValue(leftBb.psi);
                    firstFieldRight = toStringValue(rightBb.psi);
                }
            }
            if (leftBb.omega != rightBb.omega) {
                omegaDiff++;
                differs = true;
                if (firstFieldName.empty()) {
                    firstFieldName = "omega";
                    firstFieldLeft = toStringValue(leftBb.omega);
                    firstFieldRight = toStringValue(rightBb.omega);
                }
            }
            if (leftBb.n_ca_c_angle != rightBb.n_ca_c_angle) {
                nCaCDiff++;
                differs = true;
                if (firstFieldName.empty()) {
                    firstFieldName = "n_ca_c_angle";
                    firstFieldLeft = toStringValue(leftBb.n_ca_c_angle);
                    firstFieldRight = toStringValue(rightBb.n_ca_c_angle);
                }
            }
            if (leftBb.ca_c_n_angle != rightBb.ca_c_n_angle) {
                caCNDiff++;
                differs = true;
                if (firstFieldName.empty()) {
                    firstFieldName = "ca_c_n_angle";
                    firstFieldLeft = toStringValue(leftBb.ca_c_n_angle);
                    firstFieldRight = toStringValue(rightBb.ca_c_n_angle);
                }
            }
            if (leftBb.c_n_ca_angle != rightBb.c_n_ca_angle) {
                cNCaDiff++;
                differs = true;
                if (firstFieldName.empty()) {
                    firstFieldName = "c_n_ca_angle";
                    firstFieldLeft = toStringValue(leftBb.c_n_ca_angle);
                    firstFieldRight = toStringValue(rightBb.c_n_ca_angle);
                }
            }
            if (differs) {
                backboneDiffResidues++;
                if (examples.size() < 10) {
                    std::ostringstream oss;
                    oss << "residue_idx=" << j
                        << " first_field=" << firstFieldName
                        << " first_left=" << firstFieldLeft
                        << " first_right=" << firstFieldRight
                        << " res_left=" << static_cast<int>(leftBb.residue)
                        << " res_right=" << static_cast<int>(rightBb.residue)
                        << " phi=" << leftBb.phi << "/" << rightBb.phi
                        << " psi=" << leftBb.psi << "/" << rightBb.psi
                        << " omega=" << leftBb.omega << "/" << rightBb.omega
                        << " n_ca_c=" << leftBb.n_ca_c_angle << "/" << rightBb.n_ca_c_angle
                        << " ca_c_n=" << leftBb.ca_c_n_angle << "/" << rightBb.ca_c_n_angle
                        << " c_n_ca=" << leftBb.c_n_ca_angle << "/" << rightBb.c_n_ca_angle;
                    examples.push_back(oss.str());
                }
            }
        }
        std::cout << "backbone_diff_residues=" << backboneDiffResidues
                  << " residue_field_diff=" << residueDiff
                  << " phi_diff=" << phiDiff
                  << " psi_diff=" << psiDiff
                  << " omega_diff=" << omegaDiff
                  << " n_ca_c_diff=" << nCaCDiff
                  << " ca_c_n_diff=" << caCNDiff
                  << " c_n_ca_diff=" << cNCaDiff
                  << "\n";
        for (const std::string& example : examples) {
            std::cout << "backbone_diff_example " << example << "\n";
        }

        if (leftParsed.sideChainAngles != rightParsed.sideChainAngles) {
            size_t diffCount = 0;
            size_t count = std::min(leftParsed.sideChainAngles.size(), rightParsed.sideChainAngles.size());
            size_t firstDiffIndex = count;
            for (size_t j = 0; j < count; j++) {
                if (leftParsed.sideChainAngles[j] != rightParsed.sideChainAngles[j]) {
                    if (firstDiffIndex == count) {
                        firstDiffIndex = j;
                    }
                    diffCount++;
                }
            }
            std::cout << "sidechain_diff_values=" << diffCount;
            if (firstDiffIndex != count) {
                std::cout << " first_sidechain_diff_index=" << firstDiffIndex
                          << " left=" << static_cast<int>(leftParsed.sideChainAngles[firstDiffIndex])
                          << " right=" << static_cast<int>(rightParsed.sideChainAngles[firstDiffIndex]);
            }
            std::cout << "\n";
        }
        if (leftParsed.tempFactors != rightParsed.tempFactors) {
            size_t diffCount = 0;
            size_t count = std::min(leftParsed.tempFactors.size(), rightParsed.tempFactors.size());
            size_t firstDiffIndex = count;
            for (size_t j = 0; j < count; j++) {
                if (leftParsed.tempFactors[j] != rightParsed.tempFactors[j]) {
                    if (firstDiffIndex == count) {
                        firstDiffIndex = j;
                    }
                    diffCount++;
                }
            }
            std::cout << "tempfactor_diff_values=" << diffCount;
            if (firstDiffIndex != count) {
                std::cout << " first_tempfactor_diff_index=" << firstDiffIndex
                          << " left=" << static_cast<int>(leftParsed.tempFactors[firstDiffIndex])
                          << " right=" << static_cast<int>(rightParsed.tempFactors[firstDiffIndex]);
            }
            std::cout << "\n";
        }
    }

    if (leftFragments.size() != rightFragments.size()) {
        std::cout << "fragment_count_diff left=" << leftFragments.size()
                  << " right=" << rightFragments.size() << "\n";
    }
    return 0;
}

int commandContainerReport(const std::string& encodedPath) {
    std::ifstream input(encodedPath, std::ios::binary);
    if (!input) {
        std::cerr << "[Error] Failed to open " << encodedPath << "\n";
        return 1;
    }
    std::string data((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    std::string title;
    std::vector<ContainerFragment> fragments;
    bool isContainer = readContainer(data.data(), data.size(), title, fragments);
    if (!isContainer) {
        Foldcomp compRes;
        if (compRes.read(data.data(), data.size()) != 0) {
            std::cerr << "[Error] Failed to read FCZ " << encodedPath << "\n";
            return 1;
        }
        std::cout << "single_fcz"
                  << " idxResidue=" << compRes.header.idxResidue
                  << " nResidue=" << compRes.header.nResidue
                  << " idxAtom=" << compRes.header.idxAtom
                  << " nAtom=" << compRes.header.nAtom
                  << " nAnchor=" << compRes.header.nAnchor
                  << " nSideChainTorsion=" << compRes.header.nSideChainTorsion
                  << " chain=" << getChainName(compRes.header)
                  << "\n";
        std::cout << "anchor_indices";
        for (int anchorIndex : compRes.anchorIndices) {
            std::cout << " " << anchorIndex;
        }
        std::cout << "\n";
        return 0;
    }

    std::cout << "container title=" << title
              << " fragments=" << fragments.size()
              << "\n";
    for (size_t i = 0; i < fragments.size(); i++) {
        const auto& fragment = fragments[i];
        std::cout << "fragment idx=" << i
                  << " kind=" << static_cast<int>(fragment.kind)
                  << " model=" << fragment.model
                  << " chain=" << fragment.chain
                  << " payload_size=" << fragment.payload.size();
        if (fragment.kind == CONTAINER_FRAGMENT_KIND_RAW_ATOMS) {
            std::vector<AtomCoordinate> atoms;
            if (!deserializeAtomCoordinates(fragment.payload.data(), fragment.payload.size(), atoms)) {
                std::cout << " raw_decode_error=1\n";
                continue;
            }
            std::cout << " raw_atoms=" << atoms.size();
            if (!atoms.empty()) {
                std::cout << " start_resid=" << atoms.front().residue_index
                          << " end_resid=" << atoms.back().residue_index;
            }
            std::cout << "\n";
            continue;
        }
        if (fragment.kind != CONTAINER_FRAGMENT_KIND_FCZ) {
            std::cout << " unsupported_kind=" << static_cast<int>(fragment.kind) << "\n";
            continue;
        }

        Foldcomp compRes;
        if (compRes.read(fragment.payload.data(), fragment.payload.size()) != 0) {
            std::cout << " read_error=1\n";
            continue;
        }
        std::cout << " idxResidue=" << compRes.header.idxResidue
                  << " nResidue=" << compRes.header.nResidue
                  << " idxAtom=" << compRes.header.idxAtom
                  << " nAtom=" << compRes.header.nAtom
                  << " nAnchor=" << compRes.header.nAnchor
                  << " nSideChainTorsion=" << compRes.header.nSideChainTorsion
                  << " fcz_chain=" << getChainName(compRes.header)
                  << "\n";
        std::cout << "anchor_indices idx=" << i;
        for (int anchorIndex : compRes.anchorIndices) {
            std::cout << " " << anchorIndex;
        }
        std::cout << "\n";
    }
    return 0;
}

bool looksLikeDividedRoot(const std::string& sourceDir) {
    if (!fs::is_directory(sourceDir)) {
        return false;
    }
    for (const auto& entry : fs::directory_iterator(sourceDir)) {
        if (!entry.is_directory()) {
            continue;
        }
        const std::string name = entry.path().filename().string();
        if (name.size() == 2) {
            return true;
        }
    }
    return false;
}

int commandDbCompare(const std::string& sourceDir, const std::string& dbPath, bool recursive) {
    std::vector<std::string> files = getFilesInDirectory(sourceDir, recursive);
    std::sort(files.begin(), files.end());

    void* reader = make_reader(dbPath.c_str(), (dbPath + ".index").c_str(),
                               DB_READER_USE_DATA | DB_READER_USE_LOOKUP);
    if (reader == NULL) {
        std::cerr << "[Error] Failed to open database " << dbPath << "\n";
        return 1;
    }

    size_t totalFiles = 0;
    size_t proteinFiles = 0;
    size_t matchedFiles = 0;
    size_t missingEntries = 0;
    size_t decodeFailures = 0;
    size_t filesWithDiffs = 0;
    size_t totalMissingResidues = 0;
    size_t totalExtraResidues = 0;
    size_t totalChangedResidues = 0;

    for (const auto& file : files) {
        totalFiles++;
        StructureReader sourceReader;
        if (!sourceReader.load(file)) {
            continue;
        }
        std::vector<AtomCoordinate> sourceAtoms;
        sourceReader.readAllAtoms(sourceAtoms);
        removeAlternativePosition(sourceAtoms);
        if (sourceAtoms.empty()) {
            continue;
        }
        proteinFiles++;

        std::string dbKeyStr = sourcePathToDbKey(file);
        uint32_t key = reader_lookup_entry(reader, dbKeyStr.c_str());
        if (key == UINT32_MAX) {
            std::cout << "missing_db_entry " << baseName(file) << " key=" << dbKeyStr << "\n";
            missingEntries++;
            continue;
        }
        int64_t id = reader_get_id(reader, key);
        if (id == -1) {
            std::cout << "missing_db_id " << baseName(file) << " key=" << dbKeyStr << "\n";
            missingEntries++;
            continue;
        }

        std::vector<AtomCoordinate> dbAtoms;
        const char* data = reader_get_data(reader, id);
        int64_t length = reader_get_length(reader, id);
        if (data == NULL || length <= 0 || !decodeEntryToAtoms(data, static_cast<size_t>(length), dbAtoms)) {
            std::cout << "decode_failure " << baseName(file) << " key=" << dbKeyStr << "\n";
            decodeFailures++;
            continue;
        }
        removeAlternativePosition(dbAtoms);
        matchedFiles++;

        auto sourceResidues = groupResidues(sourceAtoms);
        auto dbResidues = groupResidues(dbAtoms);
        size_t missingResidues = 0;
        size_t extraResidues = 0;
        size_t changedResidues = 0;
        for (const auto& kv : sourceResidues) {
            auto it = dbResidues.find(kv.first);
            if (it == dbResidues.end()) {
                missingResidues++;
            } else if (kv.second.atoms != it->second.atoms) {
                changedResidues++;
            }
        }
        for (const auto& kv : dbResidues) {
            if (sourceResidues.find(kv.first) == sourceResidues.end()) {
                extraResidues++;
            }
        }

        totalMissingResidues += missingResidues;
        totalExtraResidues += extraResidues;
        totalChangedResidues += changedResidues;
        if (missingResidues != 0 || extraResidues != 0 || changedResidues != 0) {
            filesWithDiffs++;
            std::cout << "diff_file " << baseName(file)
                      << " missing_residue_count=" << missingResidues
                      << " extra_residue_count=" << extraResidues
                      << " changed_residue_count=" << changedResidues
                      << "\n";
        }
    }

    free_reader(reader);
    std::cout << "files_total=" << totalFiles
              << " protein_files=" << proteinFiles
              << " matched_files=" << matchedFiles
              << " missing_entries=" << missingEntries
              << " decode_failures=" << decodeFailures
              << " files_with_diffs=" << filesWithDiffs
              << " total_missing_residues=" << totalMissingResidues
              << " total_extra_residues=" << totalExtraResidues
              << " total_changed_residues=" << totalChangedResidues
              << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage:\n"
                  << "  foldcomp_audit summarize <structure>\n"
                  << "  foldcomp_audit compare <source_structure> <roundtrip_structure>\n"
                  << "  foldcomp_audit encoded-diff <left_encoded> <right_encoded>\n"
                  << "  foldcomp_audit fallback-report <structure_or_directory>\n"
                  << "  foldcomp_audit db-compare [--recursive] <source_directory> <db>\n"
                  << "  foldcomp_audit db-rmsd [--recursive] <source_directory> <db>\n"
                  << "  foldcomp_audit db-entry-compare <source_structure> <db>\n"
                  << "  foldcomp_audit dump-residues <structure>\n"
                  << "  foldcomp_audit db-entry-dump-residues <source_structure> <db>\n"
                  << "  foldcomp_audit fragment-report <structure>\n"
                  << "  foldcomp_audit trace-encode <structure> [fragment_idx]\n"
                  << "  foldcomp_audit trace-torsion-op <structure> <fragment_idx> <torsion_idx>\n"
                  << "  foldcomp_audit container-report <encoded_structure>\n";
        return 1;
    }

    std::string command = argv[1];
    if (command == "summarize") {
        return commandSummarize(argv[2]);
    }
    if (command == "compare") {
        if (argc < 4) {
            std::cerr << "[Error] compare requires two input files\n";
            return 1;
        }
        return commandCompare(argv[2], argv[3]);
    }
    if (command == "encoded-diff") {
        if (argc < 4) {
            std::cerr << "[Error] encoded-diff requires two encoded files\n";
            return 1;
        }
        return commandEncodedDiff(argv[2], argv[3]);
    }
    if (command == "fallback-report") {
        return commandFallbackReport(argv[2]);
    }
    if (command == "db-compare") {
        bool recursive = false;
        int argIndex = 2;
        if (argc >= 3 && std::string(argv[2]) == "--recursive") {
            recursive = true;
            argIndex++;
        }
        if (argc < argIndex + 2) {
            std::cerr << "[Error] db-compare requires source_directory and db path\n";
            return 1;
        }
        const std::string sourceDir = argv[argIndex];
        const std::string dbPath = argv[argIndex + 1];
        if (!recursive) {
            recursive = looksLikeDividedRoot(sourceDir);
        }
        return commandDbCompare(sourceDir, dbPath, recursive);
    }
    if (command == "db-rmsd") {
        bool recursive = false;
        int argIndex = 2;
        if (argc >= 3 && std::string(argv[2]) == "--recursive") {
            recursive = true;
            argIndex++;
        }
        if (argc < argIndex + 2) {
            std::cerr << "[Error] db-rmsd requires source_directory and db path\n";
            return 1;
        }
        const std::string sourceDir = argv[argIndex];
        const std::string dbPath = argv[argIndex + 1];
        if (!recursive) {
            recursive = looksLikeDividedRoot(sourceDir);
        }
        return commandDbRmsd(sourceDir, dbPath, recursive);
    }
    if (command == "db-entry-compare") {
        if (argc < 4) {
            std::cerr << "[Error] db-entry-compare requires source_structure and db path\n";
            return 1;
        }
        return commandDbEntryCompare(argv[2], argv[3]);
    }
    if (command == "dump-residues") {
        return commandDumpResidues(argv[2]);
    }
    if (command == "db-entry-dump-residues") {
        if (argc < 4) {
            std::cerr << "[Error] db-entry-dump-residues requires source_structure and db path\n";
            return 1;
        }
        return commandDbEntryDumpResidues(argv[2], argv[3]);
    }
    if (command == "fragment-report") {
        return commandFragmentReport(argv[2]);
    }
    if (command == "trace-encode") {
        int fragmentIndex = -1;
        if (argc >= 4) {
            fragmentIndex = std::stoi(argv[3]);
        }
        return commandTraceEncode(argv[2], fragmentIndex);
    }
    if (command == "trace-torsion-op") {
        if (argc < 5) {
            std::cerr << "[Error] trace-torsion-op requires structure, fragment_idx and torsion_idx\n";
            return 1;
        }
        return commandTraceTorsionOp(argv[2], std::stoi(argv[3]), std::stoi(argv[4]));
    }
    if (command == "container-report") {
        return commandContainerReport(argv[2]);
    }

    std::cerr << "[Error] Unknown command: " << command << "\n";
    return 1;
}
