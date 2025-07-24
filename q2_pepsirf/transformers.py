#!/usr/bin/env python
from q2_pepsirf.format_types import (
    PepsirfContingencyTSVFormat, PepsirfInfoSumOfProbesFmt,
    EnrichedPeptideDirFmt, PeptideIDListFmt, EpitopeFormat,
    MappedEpitopeFormat, GMTFormat
)
from q2_pepsirf.plugin_setup import plugin
from q2_types.feature_table import BIOMV210Format

import pandas as pd
import biom

from csv import QUOTE_NONE


# Transform a PepsirfContingencyTSVFormat into a BIOMV210Format
@plugin.register_transformer
def _0(ff: PepsirfContingencyTSVFormat) -> BIOMV210Format:
    result = BIOMV210Format()

    dataframe = pd.read_csv(str(ff), sep="\t", index_col=0)
    table = biom.Table(
        dataframe.values, observation_ids=dataframe.index,
        sample_ids=dataframe.columns
    )

    with result.open() as fh:
        table.to_hdf5(fh, generated_by="q2-pepsirf for pepsirf")

    return result

# transform a BIOMV210Format into a PepsirfContingencyTSVFormat
@plugin.register_transformer
def _1(ff: BIOMV210Format) -> PepsirfContingencyTSVFormat:
    result = PepsirfContingencyTSVFormat()

    with ff.open() as fh:
        table = biom.Table.from_hdf5(fh)
    df = table.to_dataframe(dense=True)
    df.index.name = "Sequence name"
    df.to_csv(str(result), sep="\t")

    return result

# transform a PepsirfInfoSumOfProbesFmt into a pandas dataframe
@plugin.register_transformer
def _2(ff: PepsirfInfoSumOfProbesFmt) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep="\t")
    return result

# transform a EnrichedPeptideDirFmt into a pandas dataframe
@plugin.register_transformer
def _3(ff: EnrichedPeptideDirFmt ) -> pd.DataFrame:
    pairwiseDict = {}
    for relpath, series in ff.pairwise.iter_views(pd.Series):
        pairwiseDict[relpath] = series
    df = pd.DataFrame(pairwiseDict)
    df = df.fillna(False)
    return df

# transform a PeptideIDListFmt into a pandas series
@plugin.register_transformer
def _4(ff: PeptideIDListFmt) -> pd.Series:
    with ff.open() as fh:
        ids = [id.strip() for id in fh.readlines()]
    return pd.Series(True, index = ids)

# transform a PepsirfContingencyTSVFormat into a pandas dataframe
@plugin.register_transformer
def _5(ff: PepsirfContingencyTSVFormat) -> pd.DataFrame:
    dataframe = pd.read_csv(str(ff), sep="\t", index_col=0)
    return dataframe.transpose()

# transform a PepsirfContingencyTSVFormat into a biom.Table
@plugin.register_transformer
def _6(ff: PepsirfContingencyTSVFormat) -> biom.Table:
    pass

# transform a biom.Table into a PepsirfContingencyTSV format
@plugin.register_transformer
def _7(ff: biom.Table) -> PepsirfContingencyTSVFormat:
    result = PepsirfContingencyTSVFormat()

    with open(str(result), 'w') as fh:
        ff.to_tsv(direct_io=fh, observation_column_name='Sequence name')

    return result

@plugin.register_transformer
def _9(ff: EpitopeFormat) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep='\t', index_col=0, low_memory=False)
    return result

@plugin.register_transformer
def _8(ff: pd.DataFrame) -> EpitopeFormat:
    result = EpitopeFormat()
    ff.to_csv(str(result), sep='\t')
    return result

@plugin.register_transformer
def _10(ff: MappedEpitopeFormat) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep='\t', index_col=0, low_memory=False)
    return result

@plugin.register_transformer
def _11(ff: pd.DataFrame) -> MappedEpitopeFormat:
    result = MappedEpitopeFormat()
    ff.to_csv(str(result), sep='\t')
    return result

@plugin.register_transformer
def _12(ff: GMTFormat) -> pd.DataFrame:
    pass

@plugin.register_transformer
def _13(ff: pd.DataFrame) -> GMTFormat:
    result = GMTFormat()
    ff.to_csv(str(result), sep='\t', header=False, quoting=QUOTE_NONE, escapechar=' ')
    return result
