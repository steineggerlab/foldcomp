#include "structure_codec.h"

#include "foldcomp.h"

#ifdef FOLDCOMP_WITH_STRUCTURE_READER
#include "structure_reader.h"
#endif

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <set>
#include <limits>

namespace {

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
void appendBytes(std::string& output, const T& value) {
    size_t offset = output.size();
    output.resize(offset + sizeof(T));
    memcpy(output.data() + offset, &value, sizeof(T));
}

bool isContainerHeader(const CompressedFileHeader& header) {
    bool supportedVersion = header.version == 1 || header.version == FOLDCOMP_FORMAT_VERSION_CONTAINER;
    return supportedVersion &&
           (header.flags & FOLDCOMP_FORMAT_FLAG_CONTAINER) != 0 &&
           header.nResidue == 0 &&
           header.nAtom == 0;
}

void appendTitleLines(const std::string& title, std::string& output) {
    if (title.empty()) {
        return;
    }
    const char* headerData = title.c_str();
    size_t headerLen = title.length();
    int remainingHeader = static_cast<int>(headerLen);
    char buffer[128];
    int written = snprintf(
        buffer, sizeof(buffer), "TITLE     %.*s\n",
        std::min(70, remainingHeader), headerData
    );
    if (written >= 0 && written < static_cast<int>(sizeof(buffer))) {
        output.append(buffer, written);
    }
    remainingHeader -= 70;
    int continuation = 2;
    while (remainingHeader > 0) {
        written = snprintf(
            buffer, sizeof(buffer), "TITLE  % 3d%.*s\n",
            continuation, std::min(70, remainingHeader),
            headerData + (headerLen - remainingHeader)
        );
        if (written >= 0 && written < static_cast<int>(sizeof(buffer))) {
            output.append(buffer, written);
        }
        remainingHeader -= 70;
        continuation++;
    }
}

std::string trim(const std::string& str, const std::string& whitespace = " \t") {
    const std::string::size_type strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos) {
        return "";
    }
    const std::string::size_type strEnd = str.find_last_not_of(whitespace);
    return str.substr(strBegin, strEnd - strBegin + 1);
}

bool setTitleOnFirstFczFragment(std::vector<ContainerFragment>& fragments, const std::string& title) {
    if (title.empty()) {
        return true;
    }
    for (std::vector<ContainerFragment>::iterator fragmentIt = fragments.begin();
         fragmentIt != fragments.end();
         ++fragmentIt) {
        ContainerFragment& fragment = *fragmentIt;
        if (fragment.kind != CONTAINER_FRAGMENT_KIND_FCZ) {
            continue;
        }
        const size_t minimumPayloadSize = MAGICNUMBER_LENGTH + sizeof(CompressedFileHeader);
        if (fragment.payload.size() < minimumPayloadSize) {
            return false;
        }
        if (memcmp(fragment.payload.data(), MAGICNUMBER, MAGICNUMBER_LENGTH) != 0) {
            return false;
        }

        CompressedFileHeader header = {};
        memcpy(&header, fragment.payload.data() + MAGICNUMBER_LENGTH, sizeof(header));
        const size_t anchorBytes = static_cast<size_t>(header.nAnchor) * sizeof(int);
        const size_t titleOffset = minimumPayloadSize + anchorBytes;
        if (titleOffset > fragment.payload.size()) {
            return false;
        }
        if (titleOffset + header.lenTitle > fragment.payload.size()) {
            return false;
        }

        if (title.size() > std::numeric_limits<uint32_t>::max()) {
            return false;
        }
        const uint32_t newTitleLength = static_cast<uint32_t>(title.size());
        const bool sameTitleLength = header.lenTitle == newTitleLength;
        const bool sameTitleBytes = sameTitleLength &&
            (header.lenTitle == 0 ||
             memcmp(fragment.payload.data() + titleOffset, title.data(), header.lenTitle) == 0);
        if (sameTitleBytes) {
            return true;
        }

        const size_t oldTitleEnd = titleOffset + header.lenTitle;
        header.lenTitle = newTitleLength;

        std::string updatedPayload;
        updatedPayload.reserve(fragment.payload.size() - (oldTitleEnd - titleOffset) + title.size());
        updatedPayload.append(fragment.payload.data(), MAGICNUMBER_LENGTH);
        appendBytes(updatedPayload, header);
        if (anchorBytes > 0) {
            updatedPayload.append(fragment.payload.data() + minimumPayloadSize, anchorBytes);
        }
        updatedPayload.append(title);
        updatedPayload.append(fragment.payload.data() + oldTitleEnd, fragment.payload.size() - oldTitleEnd);
        fragment.payload = std::move(updatedPayload);
        return true;
    }
    return true;
}

}

bool hasContainerMagic(const char* data, size_t size) {
    if (size < CONTAINER_MAGIC_LENGTH) {
        return false;
    }
    return memcmp(data, CONTAINER_MAGIC, CONTAINER_MAGIC_LENGTH) == 0;
}

bool parsePDBAtoms(
    const char* data, size_t size, bool singleChainOnly,
    std::vector<AtomCoordinate>& atoms, int& status
) {
    atoms.clear();
    status = PARSE_PDB_OK;
    size_t start = 0;
    std::string chain;
    while (start < size) {
        size_t end = start;
        while (end < size && data[end] != '\n') {
            end++;
        }
        size_t len = end - start;
        if (len >= 66 && memcmp(data + start, "ATOM", 4) == 0) {
            std::string line(data + start, len);
            std::string lineChain = line.substr(21, 1);
            if (singleChainOnly) {
                if (chain.empty()) {
                    chain = lineChain;
                } else if (lineChain != chain) {
                    status = PARSE_PDB_MULTIPLE_CHAINS;
                    return false;
                }
            } else {
                chain = lineChain;
            }
            try {
                atoms.emplace_back(
                    trim(line.substr(12, 4)),
                    trim(line.substr(17, 3)),
                    chain,
                    std::stoi(line.substr(6, 5)),
                    std::stoi(line.substr(22, 4)),
                    std::stof(line.substr(30, 8)),
                    std::stof(line.substr(38, 8)),
                    std::stof(line.substr(46, 8)),
                    std::stof(line.substr(54, 6)),
                    std::stof(line.substr(60, 6))
                );
            } catch (const std::exception&) {
                status = PARSE_PDB_INVALID_FORMAT;
                return false;
            }
        }
        if (end == size) {
            break;
        }
        start = end + 1;
    }
    if (atoms.empty()) {
        status = PARSE_PDB_NO_ATOM;
        return false;
    }
    return true;
}

bool parseStructureAtoms(
    const char* data, size_t size, bool singleChainOnly,
    std::vector<AtomCoordinate>& atoms, int& status, std::string* title, const char* format
) {
    atoms.clear();
    if (title != nullptr) {
        title->clear();
    }
    status = PARSE_PDB_OK;

#ifdef FOLDCOMP_WITH_STRUCTURE_READER
    StructureReader reader;
    std::string inputName;
    if (format == nullptr || format[0] == '\0' || strcmp(format, "pdb") == 0) {
        inputName = "input.pdb";
    } else if (strcmp(format, "mmcif") == 0 || strcmp(format, "cif") == 0) {
        inputName = "input.cif";
    } else {
        status = PARSE_PDB_INVALID_FORMAT;
        return false;
    }
    if (!reader.loadFromBuffer(data, size, inputName)) {
        status = PARSE_PDB_INVALID_FORMAT;
        return false;
    }

    if (title != nullptr) {
        *title = reader.title;
    }
    reader.readAllAtoms(atoms);
    if (atoms.empty()) {
        status = PARSE_PDB_NO_ATOM;
        return false;
    }

    if (singleChainOnly) {
        std::string chain;
        bool chainSet = false;
        for (const auto& atom : atoms) {
            if (!chainSet) {
                chain = atom.chain;
                chainSet = true;
            } else if (atom.chain != chain) {
                status = PARSE_PDB_MULTIPLE_CHAINS;
                atoms.clear();
                return false;
            }
        }
    }
    return true;
#else
    if (title != nullptr) {
        title->clear();
    }
    if (format != nullptr && format[0] != '\0' &&
        strcmp(format, "pdb") != 0) {
        status = PARSE_PDB_INVALID_FORMAT;
        return false;
    }
    return parsePDBAtoms(data, size, singleChainOnly, atoms, status);
#endif
}

int encodeStructureToFoldcomp(
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
    if (!parseStructureAtoms(data, size, true, atomCoordinates, status, nullptr, format)) {
        return status;
    }

    Foldcomp compRes;
    compRes.strTitle = title;
    compRes.anchorThreshold = anchorResidueThreshold;
    compRes.compress(atomCoordinates);
    if (compRes.exceedsBackboneRmsdThreshold(maxBackboneRmsd)) {
        ContainerFragment fragment;
        fragment.kind = CONTAINER_FRAGMENT_KIND_RAW_ATOMS;
        if (!atomCoordinates.empty()) {
            fragment.model = atomCoordinates.front().model;
            fragment.chain = atomCoordinates.front().chain;
        }
        if (!serializeAtomCoordinates(atomCoordinates, fragment.payload)) {
            return PARSE_PDB_INVALID_FORMAT;
        }
        std::vector<ContainerFragment> fragments;
        fragments.push_back(std::move(fragment));
        if (!writeContainerToString(output, title, fragments)) {
            return PARSE_PDB_INVALID_FORMAT;
        }
    } else {
        if (compRes.writeString(output) != 0) {
            return PARSE_PDB_INVALID_FORMAT;
        }
    }

    return PARSE_PDB_OK;
}

bool writeContainerToString(
    std::string& output, const std::string& title,
    const std::vector<ContainerFragment>& fragments
) {
    std::vector<ContainerFragment> normalizedFragments = fragments;
    if (!setTitleOnFirstFczFragment(normalizedFragments, title)) {
        return false;
    }

    output.clear();
    size_t totalSize = CONTAINER_MAGIC_LENGTH + sizeof(CompressedFileHeader) + sizeof(uint32_t) + title.size();
    for (const auto& fragment : normalizedFragments) {
        totalSize += sizeof(fragment.kind) + sizeof(fragment.model) + sizeof(uint8_t) +
                     std::min<size_t>(255, fragment.chain.size()) + sizeof(uint32_t) + fragment.payload.size();
    }
    output.reserve(totalSize);

    CompressedFileHeader header = {};
    header.version = FOLDCOMP_FORMAT_VERSION_CONTAINER;
    header.flags = FOLDCOMP_FORMAT_FLAG_CONTAINER;
    header.lenTitle = static_cast<uint32_t>(title.size());
    uint32_t nFragments = static_cast<uint32_t>(normalizedFragments.size());

    output.append(CONTAINER_MAGIC, CONTAINER_MAGIC_LENGTH);
    appendBytes(output, header);
    appendBytes(output, nFragments);
    output.append(title);

    for (const auto& fragment : normalizedFragments) {
        uint8_t chainLen = static_cast<uint8_t>(std::min<size_t>(255, fragment.chain.size()));
        uint32_t payloadSize = static_cast<uint32_t>(fragment.payload.size());
        appendBytes(output, fragment.kind);
        appendBytes(output, fragment.model);
        appendBytes(output, chainLen);
        if (chainLen > 0) {
            output.append(fragment.chain.data(), chainLen);
        }
        appendBytes(output, payloadSize);
        if (payloadSize > 0) {
            output.append(fragment.payload);
        }
    }
    return true;
}

bool readContainer(
    const char* data, size_t size, std::string& title,
    std::vector<ContainerFragment>& fragments
) {
    title.clear();
    fragments.clear();
    if (!hasContainerMagic(data, size)) {
        return false;
    }
    if (size < CONTAINER_MAGIC_LENGTH + sizeof(CompressedFileHeader) + sizeof(uint32_t)) {
        return false;
    }

    size_t offset = CONTAINER_MAGIC_LENGTH;
    CompressedFileHeader header = {};
    memcpy(&header, data + offset, sizeof(header));
    offset += sizeof(header);
    if (!isContainerHeader(header)) {
        return false;
    }

    uint32_t nFragments = 0;
    if (!readBinaryValue(data, size, offset, nFragments)) {
        return false;
    }

    if (offset + header.lenTitle > size) {
        return false;
    }
    if (header.lenTitle > 0) {
        title.assign(data + offset, data + offset + header.lenTitle);
        offset += header.lenTitle;
    }

    fragments.reserve(nFragments);
    for (uint32_t i = 0; i < nFragments; i++) {
        ContainerFragment fragment;
        fragment.kind = CONTAINER_FRAGMENT_KIND_FCZ;
        uint8_t chainLen = 0;
        uint32_t payloadSize = 0;
        if (header.version >= FOLDCOMP_FORMAT_VERSION_CONTAINER) {
            if (!readBinaryValue(data, size, offset, fragment.kind)) {
                return false;
            }
        }
        if (!readBinaryValue(data, size, offset, fragment.model) ||
            !readBinaryValue(data, size, offset, chainLen)) {
            return false;
        }
        if (offset + chainLen > size) {
            return false;
        }
        if (chainLen > 0) {
            fragment.chain.assign(data + offset, data + offset + chainLen);
            offset += chainLen;
        }
        if (!readBinaryValue(data, size, offset, payloadSize)) {
            return false;
        }
        if (offset + payloadSize > size) {
            return false;
        }
        fragment.payload.assign(data + offset, data + offset + payloadSize);
        offset += payloadSize;
        fragments.push_back(std::move(fragment));
    }
    return true;
}

void writeSegmentsToPDB(
    const std::vector<DecompressedSegment>& segments,
    const std::string& title,
    std::string& output
) {
    output.clear();
    if (segments.empty()) {
        return;
    }
    appendTitleLines(title, output);
    std::set<int> modelSet;
    for (const auto& segment : segments) {
        modelSet.insert(segment.model);
    }
    bool writeModels = modelSet.size() > 1;
    int currModel = -1;
    for (size_t i = 0; i < segments.size(); i++) {
        const auto& segment = segments[i];
        if (writeModels && segment.model != currModel) {
            if (currModel != -1) {
                output.append("ENDMDL\n");
            }
            char modelLine[32];
            int written = snprintf(modelLine, sizeof(modelLine), "MODEL     %4d\n", segment.model);
            output.append(modelLine, written);
            currModel = segment.model;
        }
        if (segment.atoms.empty()) {
            continue;
        }
        bool emitFinalTer = true;
        size_t nextIndex = i + 1;
        while (nextIndex < segments.size() && segments[nextIndex].atoms.empty()) {
            nextIndex++;
        }
        if (nextIndex < segments.size()) {
            const auto& nextSegment = segments[nextIndex];
            if (!nextSegment.atoms.empty() && nextSegment.model == segment.model) {
                const AtomCoordinate& lastAtom = segment.atoms.back();
                const AtomCoordinate& nextAtom = nextSegment.atoms.front();
                if (lastAtom.model == nextAtom.model && lastAtom.chain == nextAtom.chain) {
                    emitFinalTer = false;
                }
            }
        }
        writeAtomCoordinatesToPDB(segment.atoms, "", output, true, emitFinalTer);
    }
    if (writeModels && currModel != -1) {
        output.append("ENDMDL\n");
    }
}

#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
void writeSegmentsToMMCIF(
    const std::vector<DecompressedSegment>& segments,
    const std::string& title,
    std::string& output
) {
    std::vector<AtomCoordinate> atoms;
    size_t totalAtoms = 0;
    for (const auto& segment : segments) {
        totalAtoms += segment.atoms.size();
    }
    atoms.reserve(totalAtoms);
    for (const auto& segment : segments) {
        for (const auto& atom : segment.atoms) {
            AtomCoordinate copy = atom;
            copy.model = segment.model;
            atoms.push_back(std::move(copy));
        }
    }
    writeAtomCoordinatesToMMCIF(atoms, title, output);
}
#endif

bool decodeStructureToSegments(
    const char* data, size_t size, bool useAltOrder, std::string& title,
    std::vector<DecompressedSegment>& segments
) {
    title.clear();
    segments.clear();
    std::vector<ContainerFragment> fragments;
    if (readContainer(data, size, title, fragments)) {
        segments.reserve(fragments.size());
        for (const auto& fragment : fragments) {
            DecompressedSegment segment;
            segment.model = fragment.model;
            if (fragment.kind == CONTAINER_FRAGMENT_KIND_RAW_ATOMS) {
                if (!deserializeAtomCoordinates(fragment.payload.data(), fragment.payload.size(), segment.atoms)) {
                    return false;
                }
            } else if (fragment.kind == CONTAINER_FRAGMENT_KIND_FCZ) {
                Foldcomp compRes;
                int flag = compRes.read(fragment.payload.data(), fragment.payload.size());
                if (flag != 0) {
                    return false;
                }
                if (title.empty() && !compRes.strTitle.empty()) {
                    title = compRes.strTitle;
                }
                compRes.useAltAtomOrder = useAltOrder;
                flag = compRes.decompress(segment.atoms);
                if (flag != 0) {
                    return false;
                }
            } else {
                return false;
            }
            if (!fragment.chain.empty()) {
                for (auto& atom : segment.atoms) {
                    atom.chain = fragment.chain;
                }
            }
            segments.push_back(std::move(segment));
        }
        return true;
    }

    Foldcomp compRes;
    int flag = compRes.read(data, size);
    if (flag != 0) {
        return false;
    }
    compRes.useAltAtomOrder = useAltOrder;
    DecompressedSegment segment;
    if (compRes.decompress(segment.atoms) != 0) {
        return false;
    }
    title = compRes.strTitle;
    segments.push_back(std::move(segment));
    return true;
}

bool decodeStructureToAtoms(
    const char* data, size_t size, bool useAltOrder, std::string& title,
    std::vector<AtomCoordinate>& atoms
) {
    atoms.clear();
    std::vector<DecompressedSegment> segments;
    if (!decodeStructureToSegments(data, size, useAltOrder, title, segments)) {
        return false;
    }
    size_t totalAtoms = 0;
    for (const auto& segment : segments) {
        totalAtoms += segment.atoms.size();
    }
    atoms.reserve(totalAtoms);
    for (auto& segment : segments) {
        for (auto& atom : segment.atoms) {
            atom.model = segment.model;
            atoms.push_back(std::move(atom));
        }
    }
    return true;
}

bool decodeStructureToPDB(
    const char* data, size_t size, bool useAltOrder, std::string& title,
    std::string& pdbText
) {
    std::vector<DecompressedSegment> segments;
    if (!decodeStructureToSegments(data, size, useAltOrder, title, segments)) {
        return false;
    }
    writeSegmentsToPDB(segments, title, pdbText);
    return true;
}

#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
bool decodeStructureToMMCIF(
    const char* data, size_t size, bool useAltOrder, std::string& title,
    std::string& mmcifText
) {
    std::vector<DecompressedSegment> segments;
    if (!decodeStructureToSegments(data, size, useAltOrder, title, segments)) {
        return false;
    }
    writeSegmentsToMMCIF(segments, title, mmcifText);
    return true;
}
#endif
