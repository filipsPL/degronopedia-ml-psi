#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import peptides
from rdkit.Chem import AllChem
# from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import DataStructs

import seaborn as sns

import numpy as np
import pandas as pd
import argparse
import os

from catboost import CatBoostRegressor

# sequence="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRT"

columnSet = "^Whole_.*|^Ogryzek4_.*‎"

nBits = 128

descriptors = [
    'qed', 'MolWt', 'MaxPartialCharge', 'MinPartialCharge', 'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3',
    'BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v',
    'Chi4n', 'Chi4v', 'HallKierAlpha', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'PEOE_VSA1', 'PEOE_VSA10',
    'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5',
    'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4',
    'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11',
    'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8',
    'SlogP_VSA9', 'TPSA', 'EState_VSA1', 'EState_VSA10', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4',
    'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'VSA_EState1', 'VSA_EState10',
    'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8',
    'VSA_EState9', 'FractionCSP3', 'MolLogP', 'MolMR'
]
fpHeaders = ["morganFp2_%i" % (i) for i in range(0, nBits)]

# --------------- #

allHeaders = []
allHeaders.extend(descriptors)
allHeaders.extend(fpHeaders)

# --------------- #

calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors)

# --------------- #


# select subset of columns from the big set of columns, basing on the regex
def selectColumnSubset(df, columnSet):
    outDf = df.filter(regex=(columnSet))
    return outDf


def file_read(filename):
    f = open(filename, 'r')
    content = f.read()
    f.close()
    return content


def generate_ogryzki(sequence):

    ogryzki_C = []
    ogryzki_N = []
    if sequence[0] == 'M':
        sequence = sequence[1:]

    # first/last 2, 3 ... 10 residues
    for i in range(2, 11):
        s_c = sequence[-i:]
        s_n = sequence[:i]
        ogryzki_C.append(s_c)
        ogryzki_N.append(s_n)

    return [ogryzki_C, ogryzki_N]


def generate_peptides(sequence, ogryzki):

    mega_dict = {}
    mega_dict['query'] = []
    window_columns = []

    whole = peptides.Peptide(sequence).descriptors()
    for key in whole.keys():
        mega_dict['query'].append(whole[key])
        window_columns.append('Whole_peptides_' + key)

    for og in range(len(ogryzki)):
        x = peptides.Peptide(ogryzki[og]).descriptors()
        for kx in x.keys():
            mega_dict['query'].append(x[kx])
            window_columns.append(f'Ogryzek{og+2}_peptides_' + kx)

    df = pd.DataFrame.from_dict(mega_dict, orient='index', columns=window_columns)
    df.index.name = 'Index'

    return df


def aa_gravy(sequence, ogryzki):

    mega_dict2 = {}
    mega_dict2['query'] = []

    for i in range(len(sequence)):
        mega_dict2['query'].append(sequence[i])
    whole = round(ProteinAnalysis(sequence).gravy(), 2)
    mega_dict2['query'].append(whole)

    for og in range(len(ogryzki)):
        for i in range(len(ogryzki[og])):
            mega_dict2['query'].append(ogryzki[og][i])
        x = round(ProteinAnalysis(ogryzki[og]).gravy(), 2)
        mega_dict2['query'].append(x)

    cols = []
    for b in range(len(sequence)):
        cols.append(f'Whole_aa_{b+1}')
    cols.append('Whole_gravy')

    for i in range(2, 11):
        for b in range(i):
            cols.append(f'Ogryzek{i}_aa_{b+1}')
        cols.append(f'Ogryzek{i}_gravy')

    df_aa = pd.DataFrame.from_dict(mega_dict2, orient='index', columns=cols)
    df_aa.index.name = 'Index'

    return df_aa


def computeMorganFP(mol, depth=2, nBits=nBits):
    a = np.zeros(nBits, dtype=int)
    try:
        DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(mol, depth, nBits), a)
    except:
        return None
    return a


def calcRDkitDescs(m):
    ds = list(calc.CalcDescriptors(m))
    return ds


def calcDescForSeq(seq):

    m = AllChem.rdmolfiles.MolFromFASTA(seq)

    ## output
    outVector = []

    ## descriptors
    ds = calcRDkitDescs(m)
    outVector.extend(ds)

    ## fingerprint
    fp = computeMorganFP(m)
    outVector.extend(fp)

    return outVector


def generate_rdkit(sequence, ogryzki):
    filip_mega_dict = {}
    filip_mega_dict['query'] = []

    whole = calcDescForSeq(sequence)

    filip_mega_dict['query'].extend(whole)

    for og in range(len(ogryzki)):
        x = calcDescForSeq(ogryzki[og])
        filip_mega_dict['query'].extend(x)

    ### columns ###

    filip_columns = []
    for c in allHeaders:
        filip_columns.append('Whole_RDKit_' + c)

    for i in range(2, 11):
        for c in allHeaders:
            filip_columns.append(f'Ogryzek{i}_RDKit_' + c)

    df_filip = pd.DataFrame.from_dict(filip_mega_dict, orient='index', columns=filip_columns)
    df_filip.index.name = 'Index'

    return df_filip


def generate_ML_tsv(seq, ogryzki):

    df1 = generate_peptides(seq, ogryzki)
    df2 = aa_gravy(seq, ogryzki)
    df3 = generate_rdkit(seq, ogryzki)

    df_aa_ML = df1.merge(df2, on='Index')
    df_aa_ML = df_aa_ML.merge(df3, on='Index')
    # df_aa_ML.to_csv(filename + '.tsv', sep='\t')

    return df_aa_ML


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--sequence', dest='inputFile', required=True, help='file with the sequence in plain text')

    parser.add_argument('--type',
                        dest='type',
                        required=True,
                        choices=['C', 'NiMetNo', 'NiMetYes'],
                        help='prediction type to make')

    args = parser.parse_args()
    sequenceFile = args.inputFile
    type = args.type

    if type == 'C':
        print("C-terminus")
        modelFile = "c-terminus.cbm"
        ogryzkiListNo = 0

    elif type == 'NiMetNo':
        print("N-terminus with initiator Met cleaved")
        modelFile = "n-terminus_iMetNO.cbm"
        ogryzkiListNo = 1

    elif type == 'NiMetYes':
        print("N-terminus with initiator Met NOT cleaved")
        modelFile = "n-terminus_iMetYES.cbm"
        ogryzkiListNo = 1

    else:
        print("Mission imposible")
        exit(2)

    # read sequence here
    f = open(sequenceFile, 'r')
    sequence = f.read().strip()
    f.close()

    if type == 'C':
        seq2 = sequence[-23:]
    elif sequence[0] == 'M':
        print("The sequence contains M at the first position")
        seq2 = sequence[1:24]
    else:
        seq2 = sequence[:23]

    try:
        ogryzkiList = generate_ogryzki(
            sequence)  # currently, iMetNO/YES yields exactly the same descriptors - their values are identical

    except Exception:
        ogryzkiList = [None, None]

    # input data with descriptors
    X = generate_ML_tsv(seq2, ogryzkiList[ogryzkiListNo])

    # set up the regressor
    reg = CatBoostRegressor()

    # read model file regardless of the working directory
    model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models", modelFile)
    reg.load_model(model_path, format='cbm')

    # preds
    X = X[reg.feature_names_]
    preds = reg.predict(X)

    print("Predicted PSI: %.2f" % (preds[0]))

    exit(0)
