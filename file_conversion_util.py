import base64
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor,AllChem
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import rdMolDraw2D
from itertools import chain
from io import StringIO, BytesIO
import re


SDF_FILE = 'SDF_FILE'
CSV_FILE = 'CSV_FILE'
XLS_FILE = 'XLS_FILE'
TXT_FILE = 'TXT_FILE'

STRUCTURE = 'Molecule'
SMILES = 'SMILES'

DO_NOT_HIGHLIGHT = "Do not highlight"
HIGHLIGHT_ALL = "Highlight All"
HIGHLIGHT_UNIQUE = "Highlight Unique"



def get_file_type(uploaded_file):
    if uploaded_file is not None:
        if uploaded_file.name.upper().endswith('.CSV'):
            return CSV_FILE
        elif uploaded_file.name.upper().endswith('.XLS') or uploaded_file.name.upper().endswith('.XLSX'):
            return XLS_FILE
        elif uploaded_file.name.upper().endswith('.SDF') or uploaded_file.name.upper().endswith('.SD'):
            return SDF_FILE
        elif uploaded_file.name.upper().endswith('.TXT'):
            return TXT_FILE
        else:
            return TXT_FILE
    else:
        return None


def get_sdf_df(uploaded_file):
    if uploaded_file is not None:
        df_upload = PandasTools.LoadSDF(uploaded_file, idName='__ID__', smilesName=SMILES, molColName=STRUCTURE)
        del df_upload['__ID__']
        return df_upload
    

def get_list(inputs):
    input_list = []
    if inputs:
        input_list = re.split(',|\n', inputs)
        input_list = [input for input in input_list if input.strip()]
    return input_list

def get_rename_dict(header_rename):
    if header_rename:
        rename_list = re.split(',|\n', header_rename)
        rename_list = [rename for rename in rename_list if rename.strip()]
        #rename_list = header_rename.split(',')
        rename_dict = {}
        for rename in rename_list:
            rename_old_new = rename.split('=')
            if len(rename_old_new) ==2:
                rename_dict[rename_old_new[0].strip()] = rename_old_new[1].strip()

        return rename_dict

# def get_csv_headers(upload_file):
#     f = StringIO()
#
#     unique_headers = set()
#     with open(upload_file, 'rb') as fin:
#         csvin = csv.reader(fin)
#         unique_headers.update(next(csvin, []))

def get_value_from_key(key, kv_dict):
    if key and key in kv_dict:
        return kv_dict.get(key)
    else:
        return key

def get_df_xlsx(df, sheet_name=None):
    f = BytesIO()
    if sheet_name:
        df.to_excel(f, sheet_name=sheet_name, index=False)
    else:
        df.to_excel(f, index=False)
        
    return f.getvalue()

def get_df_csv(df, sep):
    f = StringIO()
    df.to_csv(f, sep=sep, index=False)
    return f.getvalue()

def get_df_sdf(df, structure_column):
    f = StringIO()
    PandasTools.WriteSDF(df, f, molColName=structure_column, properties=list(df.columns), allNumeric=False)
    return f.getvalue()

def get_val_list(df, col):
    all_list = df[col].dropna().tolist()
    distinct_list = list(set(all_list))
    distinct_list.sort()
    distinct_list.insert(0, '--')
    return distinct_list


def split_dataframe(df, chunk_size):
    chunks = list()
    num_chunks = (len(df)-1)//chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
    return chunks


def moltosvg(mol, molSize = (800,400), kekulize = False, highlight_sub=None, highlight_mode=DO_NOT_HIGHLIGHT):
    
    if  highlight_sub == None: # Cannot highlight if highlight_sub not provided
        highlight_mode=DO_NOT_HIGHLIGHT
    
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    
    if highlight_mode == DO_NOT_HIGHLIGHT:
        drawer.DrawMolecule(mc)
    elif highlight_mode in (HIGHLIGHT_UNIQUE, HIGHLIGHT_ALL):
        highlight_tt = mc.GetSubstructMatches(highlight_sub)
        hightlight_shape = np.shape(highlight_tt)
        if hightlight_shape[0] == 1:
            highlight_tuple = tuple(chain.from_iterable(highlight_tt))
            drawer.DrawMolecule(mc, highlightAtoms=highlight_tuple)
        else:
            if highlight_mode == HIGHLIGHT_UNIQUE:
                drawer.DrawMolecule(mc)
            elif highlight_mode == HIGHLIGHT_ALL:
                highlight_tuple = tuple(chain.from_iterable(highlight_tt))
                drawer.DrawMolecule(mc, highlightAtoms=highlight_tuple)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg.replace('svg:','')
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = rf'<img src="data:image/svg+xml;base64, {b64}"/>'
    return html