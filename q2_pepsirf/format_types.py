#!/usr/bin/env python
from q2_types.feature_table import FeatureTable
from qiime2.plugin import SemanticType

import os
import qiime2.plugin.model as model

# create a semantic type for each format type created
Normed = SemanticType("Normed", variant_of=FeatureTable.field["content"])
NormedDifference = SemanticType(
    "NormedDifference",
    variant_of=FeatureTable.field["content"]
)
NormedDiffRatio = SemanticType(
    "NormedDiffRatio",
    variant_of=FeatureTable.field["content"]
)
NormedRatio = SemanticType(
    "NormedRatio",
    variant_of=FeatureTable.field["content"]
)
NormedSized = SemanticType(
    "NormedSized",
    variant_of=FeatureTable.field["content"]
)
Epitope = SemanticType(
    "Epitope",
    variant_of=FeatureTable.field["content"]
)
MappedEpitope = SemanticType(
    "MappedEpitope",
    variant_of=FeatureTable.field["content"]
)
Zscore = SemanticType("Zscore", variant_of=FeatureTable.field["content"])
RawCounts = SemanticType("RawCounts", variant_of=FeatureTable.field["content"])
PairwiseEnrichment = SemanticType("PairwiseEnrichment")
ZscoreNan = SemanticType("ZscoreNan")
PeptideBins = SemanticType("PeptideBins")
InfoSNPN = SemanticType("InfoSNPN")
InfoSumOfProbes = SemanticType("InfoSumOfProbes")
EnrichThresh = SemanticType("EnrichThresh")
MultiFile = SemanticType("MultiFile")
ProteinFasta = SemanticType("ProteinFasta")
PeptideFasta = SemanticType("PeptideFasta")
Link = SemanticType("Link")
PeptideIDList = SemanticType("PeptideIDList")
PepsirfDMP = SemanticType("PepsirfDMP")
DeconvSingluar = SemanticType("DeconvSingluar")
ScorePerRound = SemanticType("ScorePerRound")
DeconvBatch = SemanticType("DeconvBatch")
PeptideAssignmentMap = SemanticType("PeptideAssignmentMap")
DemuxFif = SemanticType("DemuxFif")
DemuxSampleList = SemanticType("DemuxSampleList")
DemuxIndex = SemanticType("DemuxIndex")
DemuxLibrary = SemanticType("DemuxLibrary")
DemuxFastq = SemanticType("DemuxFastq")
DemuxDiagnostic = SemanticType("DemuxDiagnostic")
ProteinAlignment = SemanticType("ProteinAlignment")
MutantReference = SemanticType("MutantReference")

# create a format for epitope metadata
class EpitopeFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        pass

EpitopeDirFmt = model.SingleFileDirectoryFormat(
    "EpitopeDirFmt", "epitope.tsv", EpitopeFormat
)

class MappedEpitopeFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        pass

MappedEpitopeDirFmt = model.SingleFileDirectoryFormat(
    "MappedEpitopeDirFmt", "mapped-epitope.tsv", MappedEpitopeFormat
)

# create a format for a featuretable file
class PepsirfContingencyTSVFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            # Read lines skipping comments until we either hit a valid header
            # line and return or hit an invalid one and fail out.
            while True:
                line = fh.readline()

                # Skip comments
                if line.startswith('#'):
                    continue

                # The first valid line should start with this
                if not line.startswith("Sequence name\t"):
                    raise model.ValidationError(
                        'TSV does not start with "Sequence name"')

                return


PepsirfContingencyTSVDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfContingencyTSVDirFmt", "pepsirf-table.tsv",
    PepsirfContingencyTSVFormat
)


# create a format for a peptide id list file
class PeptideIDListFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == []:
                raise model.ValidationError(
                    "TXT file is empty, not a enriched peptides file."
                )


# create a format for enrichment failure
class EnrichmentFailureFmt(model.TextFileFormat):
    def _validate_( self, level="min"):
        pass


# create a format for enrichment directory
class EnrichedPeptideDirFmt(model.DirectoryFormat):
    pairwise = model.FileCollection(
        r".+_+.+\.txt",
        format=PeptideIDListFmt
    )
    @pairwise.set_path_maker
    def pairwise_pathmaker(self, comparisons, suffix):
        return f'{"~".join(comparisons)}_{suffix}.txt'

    failures = model.File(
        "failedEnrichment.txt",
        format=EnrichmentFailureFmt
    )


PeptideIDListDirFmt = model.SingleFileDirectoryFormat(
    "PeptideIDListDirFmt", "samp_A~samp_B.txt",
    PeptideIDListFmt
)


# create a format for zscore nan file
class ZscoreNanFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("Probe name\t"):
                    raise model.ValidationError(
                        'nan does not contain "Probe name" as the first'
                        " header. Not a .nan file."
                    )


ZscoreNanDirFmt = model.SingleFileDirectoryFormat(
    "ZscoreNanDirFmt", "nan-zscores.nan",
    ZscoreNanFormat
)


# create a format for bins file
class PeptideBinFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError("TSV is empty, not a bins file.")


PeptideBinDirFmt = model.SingleFileDirectoryFormat(
    "PeptideBinDirFmt", "bins.tsv",
    PeptideBinFormat
)


# create a format for info num of samples, num of probes files
class PepsirfInfoSNPNFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError("TSV is empty.")


PepsirfInfoSNPNDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfInfoSNPNDirFmt", "SN.tsv",
    PepsirfInfoSNPNFormat
)


# create a file for info Sum of probes file
class PepsirfInfoSumOfProbesFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("Sample name\t"):
                    raise model.ValidationError(
                        'TSV does not start with "Sample name"'
                    )


PepsirfInfoSumOfProbesDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfInfoSumOfProbesDirFmt", "RC.tsv",
    PepsirfInfoSumOfProbesFmt
)


# create a format for a threshold file
class EnrichThreshFileFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError("TSV is empty.")


EnrichThreshFileDirFmt = model.SingleFileDirectoryFormat(
    "EnrichThreshFileDirFmt", "_thresh.tsv",
    EnrichThreshFileFormat
)


# create a format for a multifile
class SubjoinMultiFileFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError(".txt file is empty.")


SubjoinMultiFileDirFmt = model.SingleFileDirectoryFormat(
    "SubjoinMultiFileDirFmt", "multifile.txt",
    SubjoinMultiFileFmt
)


# create a format for a protein fasta file
class ProteinFastaFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError(".fasta file is empty.")


ProteinFastaDirFmt = model.SingleFileDirectoryFormat(
    "ProteinFastaDirFmt", "protein_file.fasta",
    ProteinFastaFmt
)


# create a format for a peptide fasta file
class PeptideFastaFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError(".faa file is empty.")


PeptideFastaDirFmt = model.SingleFileDirectoryFormat(
    "PeptideFastaDirFmt", "peptide_file.faa",
    PeptideFastaFmt
)


# create a format for a link file
class PepsirfLinkTSVFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("Peptide Name\t"):
                    raise model.ValidationError(
                        'TSV does not start with "Peptide Name"'
                    )


PepsirfLinkTSVDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfLinkTSVDirFmt", "link_output.tsv",
    PepsirfLinkTSVFormat
)


# create a format for a DMP file
class PepsirfDMPFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError("TSV is empty, not a .dmp file.")


PepsirfDMPDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfDMPDirFmt", "file.dmp",
    PepsirfDMPFormat
)


# create a format for a deoncv singular file
class PepsirfDeconvSingularFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("Species Name\t"):
                    raise model.ValidationError(
                        'TSV does not start with "Species Name"'
                    )


PepsirfDeconvSingularDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfDeconvSingularDirFmt", "singular_mode.tsv",
    PepsirfDeconvSingularFormat
)


# create a format for a deconv batch dir
class PepsirfDeconvBatchDirFmt(model.DirectoryFormat):
    batch = model.FileCollection(
        r".+_+.+\.txt",
        format=PepsirfDeconvSingularFormat
    )
    @batch.set_path_maker
    def batch_pathmaker(self, comparisons, suffix):
        return f'{"~".join(comparisons)}_{suffix}.txt'


# create a format for a score-per-round dir
class ScorePerRoundFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("Species ID\t"):
                    raise model.ValidationError(
                        'round file does not start with "Species ID"'
                    )


class ScorePerRoundDirFmt(model.DirectoryFormat):
    score = model.FileCollection(r"round_.+", format=ScorePerRoundFmt)
    @score.set_path_maker
    def score_pathmaker(self, round):
        return f"round_{round}"


# create a format for a peptide assignment map dir
class PeptideAssignMapFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("Peptide\t"):
                    raise model.ValidationError(
                        '.map does not start with "Peptide"'
                    )


class PeptideAssignMapDirFmt(model.DirectoryFormat):
    batch = model.FileCollection(r".+_.+\.map", format=PeptideAssignMapFormat)
    @batch.set_path_maker
    def batch_pathmaker(self, comparisons, suffix):
        return f'{"~".join(comparisons)}_{suffix}.map'


# create a format for a demux fif file
class PepsirfDemuxFifFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError(
                    "TSV is empty, not a demux fif file."
                )


PepsirfDemuxFifDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfDemuxFifDirFmt", "demux_fif.tsv",
    PepsirfDemuxFifFmt
)


# create a format for a demux sample list file
class PepsirfDemuxSampleListFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError(
                    "TSV is empty, not a demux sample list file."
                )


PepsirfDemuxSampleListDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfDemuxSampleListDirFmt", "demux_index.tsv",
    PepsirfDemuxSampleListFmt
)


# create a format for a demux index file
class PepsirfDemuxIndexFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError(
                    ".fa is empty, not a demux index file."
                )


PepsirfDemuxIndexDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfDemuxIndexDirFmt", "demux_index.fa",
    PepsirfDemuxIndexFmt
)


# create a format for a demux library file
class PepsirfDemuxLibraryFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError(
                    ".fna is empty, not a demux library file."
                )


PepsirfDemuxLibraryDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfDemuxLibraryDirFmt", "demux_library.fna",
    PepsirfDemuxLibraryFmt
)


# create a format for a demux fastq file
class PepsirfDemuxFastqFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            line = list(zip(range(1), fh))
            if line == [] :
                raise model.ValidationError(
                    ".fna is empty, not a demux library file."
                )


PepsirfDemuxFastqDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfDemuxFastqDirFmt", "demux.fastq",
    PepsirfDemuxFastqFmt
)


# create a format for a demux diagnostic file
class PepsirfDemuxDiagnosticFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("Sample name\t"):
                    raise model.ValidationError(
                        'TSV does not start with "Sample name"'
                    )


PepsirfDemuxDiagnosticDirFmt = model.SingleFileDirectoryFormat(
    "PepsirfDemuxDiagnosticDirFmt", "diagnostic_output.tsv",
    PepsirfDemuxDiagnosticFormat
)


# create a format for a protein alignment file
class ProteinAlignmentManifestFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("ProtName\t"):
                    raise model.ValidationError(
                        'TSV does not start with "ProtName"'
                    )


class PeptideToProteinAlignmentFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        #validate header
        pass


# class ProteinAlignmentDirFormat(model.DirectoryFormat):
#     manifest=model.File(
#         "manifest.tsv",
#         format=ProteinAlignmentManifestFormat
#     )
#     alignments=model.FileCollection(
#         r".+Aligned\.txt",
#         format=PeptideToProteinAlignmentFormat,
#     )
#     @alignments.set_path_maker
#     def alignment_pathmaker(self, name):
#         return "%sAligned.txt" % name


class ProteinAlignmentFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        pass


ProteinAlignmentDirFormat = model.SingleFileDirectoryFormat(
    "ProteinAlignmentDirFormat", "manifest.tsv",
    ProteinAlignmentFmt
)


# create a format for a demux diagnostic file
class MutantReferenceFileFmt(model.TextFileFormat):
    def _validate_(self, level="min"):
        with self.open() as fh:
            for _, line in zip(range(1), fh):
                if not line.startswith("CodeName\t"):
                    raise model.ValidationError(
                        'TSV does not start with "CodeName"'
                    )


MutantReferenceDirFmt = model.SingleFileDirectoryFormat(
    "MutantReferenceDirFmt", "references.tsv",
    MutantReferenceFileFmt
)

