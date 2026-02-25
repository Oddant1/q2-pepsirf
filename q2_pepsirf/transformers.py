#!/usr/bin/env python
import ast

from q2_pepsirf.format_types import (
    PepsirfContingencyTSVFormat, PepsirfInfoSumOfProbesFmt,
    EnrichedPeptideDirFmt, PeptideIDListFmt, EpitopeFormat,
    MappedEpitopeFormat, GMTFormat, EnrichedFormat, PSEAScoresFormat
)
from q2_pepsirf.plugin_setup import plugin
from q2_types.feature_table import BIOMV210Format

import pandas as pd
import biom


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

@plugin.register_transformer
def _6(ff: EpitopeFormat) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep='\t', index_col=0, low_memory=False)
    return result

@plugin.register_transformer
def _7(ff: pd.DataFrame) -> EpitopeFormat:
    result = EpitopeFormat()
    ff.to_csv(str(result), sep='\t')
    return result

@plugin.register_transformer
def _8(ff: MappedEpitopeFormat) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep='\t', index_col=0, low_memory=False)
    result['Subtype'] = result['Subtype'].apply(ast.literal_eval)
    result['CodeName'] = result['CodeName'].apply(ast.literal_eval)
    return result

@plugin.register_transformer
def _9(ff: pd.DataFrame) -> MappedEpitopeFormat:
    result = MappedEpitopeFormat()
    ff.to_csv(str(result), sep='\t')
    return result

@plugin.register_transformer
def _10(ff: GMTFormat) -> pd.DataFrame:
     result = pd.DataFrame(columns=['EpitopeID'])

     with open(str(ff)) as fh:
         for line in fh.readlines():
             speciesID, epitopeID = line.split('\t\t')
             epitopeID = epitopeID.split('\t')
             result.loc[speciesID] = [epitopeID]

     result.index.name = 'SpeciesID'
     return result

@plugin.register_transformer
def _11(ff: pd.DataFrame) -> GMTFormat:
    result = GMTFormat()

    with open(str(result), 'w') as fh:
        for _, row in ff.iterrows():
            line = row.name + "\t\t"

            for elem in row['EpitopeID']:
                line += elem + '\t'

            line = line.rstrip()
            line += '\n'
            fh.write(line)

    return result

@plugin.register_transformer
def _12(ff: pd.DataFrame) -> EnrichedFormat:
    result = EnrichedFormat()
    ff.to_csv(str(result), sep='\t')
    return result

@plugin.register_transformer
def _13(ff: pd.DataFrame) -> PSEAScoresFormat:
    result = PSEAScoresFormat()
    ff.to_csv(str(result), sep='\t')
    return result

@plugin.register_transformer
def _14(ff: PSEAScoresFormat) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep='\t', index_col=0, low_memory=False)
    return result
