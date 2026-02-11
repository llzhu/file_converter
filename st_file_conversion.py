import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
import os.path
from os import path
import streamlit as st
from datetime import date
import math
from file_conversion_util import *




st.set_page_config(page_title='Data file conversion', layout='wide')


 
st.markdown("**This is a tool for processing and converting csv, excel and sd files**")
# st.badge("Structure can be displayed for sd file, and csv/excel file if SMILES field exists.")

header_expander = st.expander('Download')
with header_expander:
    c1, c2, c3 = st.columns([3, 3, 3])
    with c1:
        excel_download_container = st.container()
        # mapping_container = st.container()
    with c2:
        csv_download_container = st.container()
        # merge_container = st.container()
    with c3:
        sdf_download_container = st.container()
        
    # with c4:
    #     headers_container = st.container()



details_expander = st.sidebar.expander('Structure:', expanded=True)

upload_expander = st.sidebar.expander("Upload Input File")
df_upload = None
upload_file_name = ''
delimiter = ','  # default delimiter for csv file
with upload_expander:
    uploaded_file = upload_expander.file_uploader("Input file:")
    if uploaded_file:
        upload_file_name = uploaded_file.name
        upload_file_name = os.path.splitext(upload_file_name)[0]
       
        file_type = get_file_type(uploaded_file)
        
        if file_type == CSV_FILE or file_type == TXT_FILE:
            delimiter = st.selectbox('Delimitor:', [',', 'tab', ';', 'space'])
            if delimiter == 'tab':
                delimiter = '\t'
            elif delimiter == 'space':
                delimiter = ' '

        elif file_type == XLS_FILE:
            xls_f = pd.ExcelFile(uploaded_file)
            sheet_names = xls_f.sheet_names
            if len(sheet_names) > 1:
                sheet_names.insert(0, '--')
                sheet_name = st.selectbox('Sheet name:', sheet_names, key='mapping_sheet')
            else:
                sheet_name = 0

            skip_rows = st.number_input('Skip Rows', min_value=0, max_value=10, key='skip')

        if file_type in (CSV_FILE, TXT_FILE) and delimiter != '--':
            forced_header_str = st.text_area('Forced string columns: (separate by , or newline')
            forced_header_list = get_list(forced_header_str)

            dict_dtypes = None
            if forced_header_list and len(forced_header_list) > 0:
                dict_dtypes = {x: str for x in forced_header_list}
            df_upload = pd.read_csv(uploaded_file, sep=delimiter, dtype=dict_dtypes)

        elif file_type == XLS_FILE and sheet_name != '--':
            forced_header_str = st.text_area('Forced string columns: (separate by , or newline')
            forced_header_list = get_list(forced_header_str)

            dict_dtypes = None
            if forced_header_list and len(forced_header_list) > 0:
                dict_dtypes = {x: str for x in forced_header_list}
            df_upload = pd.read_excel(uploaded_file, sheet_name=sheet_name, skiprows=skip_rows, dtype=dict_dtypes, usecols=lambda x: 'Unnamed' not in x)
            df_upload = df_upload.dropna(how='all')
        elif file_type == SDF_FILE:
            df_upload = get_sdf_df(uploaded_file)
        else:
            st.stop()
    else:
        st.stop()

if df_upload is not None:
    headers = list(df_upload)
else:
    st.stop()

#
preprocess_expander = st.sidebar.expander("Pre-processing")
with preprocess_expander:
    headers = st.multiselect('Select headers', headers, headers, key='headers_pre')

    header_replace = st.text_area('Relace part of Header names: (old1=new1,old2=new2, ....)', key='replace')
    replace_dict = get_rename_dict(header_replace)

    if replace_dict and len(replace_dict) > 0:
        df_upload = df_upload[headers]
        new_headers = headers.copy()
        for key, value in replace_dict.items():
            new_headers = [h.replace(key, value) for h in new_headers]

        df_upload = df_upload.rename(columns=dict(zip(headers, new_headers)))
    else:
        df_upload = df_upload[headers]

    headers = list(df_upload)
    headers_select = headers.copy()
    headers_select.insert(0, '--')

    num_merge = st.number_input('Numbers of pairs to merge:', value=0, key='num_merge')
    if num_merge > 0:
        for group in range(num_merge):
            header_1 = st.selectbox(f'Select header 1 for group {group+1}', headers_select, key=f'header1{group}')
            header_2 = st.selectbox(f'Select header 2 for group {group+1}', headers_select, key=f'header2{group}')
            sep = st.text_input('Separated by', value='-', key=f'new{group}')
            header_new = st.text_input('New header name', value='', key=f'new{group}')
            if header_1 != '--' and header_2 != '--' and header_new:
                df_upload[header_new] = df_upload[header_1].astype(str) + sep + df_upload[header_2].astype(str)

                del df_upload[header_1]
                del df_upload[header_2]

                headers = list(df_upload)
                headers_select = headers.copy()
                headers_select.insert(0, '--')


pivot_expander = st.sidebar.expander("Pivot")
with pivot_expander:
    all_headers = list(df_upload)
    id_headers = st.multiselect('Select id headers', all_headers, key='p_id_headers')
    pivot_on = st.multiselect('Select headers to pivot on', all_headers, key='p_pivot_on')
    value_header = st.selectbox('Select value header', all_headers, key='p_value_header')
    
    pivot_go = st.checkbox('Go', key='pivot_go')

    if id_headers and pivot_on and pivot_go:
        df_upload = df_upload.pivot_table(index=id_headers, columns=pivot_on, values=value_header, aggfunc='first')
        df_upload = df_upload.reset_index() 

headers = list(df_upload)



unpivot_expander = st.sidebar.expander("Unpivot")
with unpivot_expander:
    all_headers = list(df_upload)
    id_headers = st.multiselect('Select id headers', all_headers, key='u_id_headers')
    melt_headers = st.multiselect('Select unpivot headers', all_headers, key='u_melt_headers')
    
    unpivot_go = st.checkbox('Go', key='unpivot_go')

    if id_headers and melt_headers and unpivot_go:
        df_upload = pd.melt(df_upload, id_vars=id_headers, value_vars=melt_headers,
            var_name='value_desc', value_name='value')
        df_upload = df_upload[df_upload['value'].notnull()]
        # df_upload = df_upload[df_upload['value'].str.len() > 0]
headers = list(df_upload)


add_field_expander = st.sidebar.expander('Create new fields by parsing another field')
with add_field_expander:
    num_parse = st.number_input('Numbers of fields to parse:', value=0, key='nf_parse')
    if num_parse:
        for i_parse in range(num_parse):
            parsing_from = st.selectbox(f'Parse from column {i_parse + 1}:', list(df_upload), key=f'field_{i_parse}')
            parsing_re = st.text_input(f"Parsing Reg Expression for column {i_parse + 1}:", key=f're_{i_parse}')
            pattern = re.compile(parsing_re, re.IGNORECASE)
            n_group = pattern.groups
            field_name_re = ['' for i in range(n_group + 1)]
            if parsing_from and parsing_re and n_group > 0:
                for g in range(1,n_group+1):
                    field_name_re[g] = st.text_input(f'Column name for group {g}', key=f'fn_re_{i_parse}_{g}')
                    if field_name_re[g]:
                        df_upload[field_name_re[g]] = df_upload[parsing_from].apply(lambda x: pattern.match(x).group(g))


headers = list(df_upload)

merge_expander = st.sidebar.expander("Merge External File")
with merge_expander:
    df_merge = None
    merge_file = st.file_uploader("File to be merged:")
    if merge_file:
        file_type = get_file_type(merge_file)
        if file_type == CSV_FILE:
            delimiter = st.sidebar.selectbox('Delimitor:', [',', 'tab', ';'], key='merge_del')
            if delimiter == 'tab':
                delimiter = '\t'
        elif file_type == XLS_FILE:
            xls_f = pd.ExcelFile(merge_file)
            sheet_names = xls_f.sheet_names
            if len(sheet_names) > 1:
                sheet_names.insert(0, '--')
                sheet_name = st.selectbox('Sheet name:', sheet_names, key='merge_sheet')
            else:
                sheet_name = 0
            skip_rows = st.number_input('Skip Rows', min_value=0, max_value=10, key='merge_skip')

        if file_type == CSV_FILE and delimiter != '--':
            df_merge = pd.read_csv(merge_file, sep=delimiter)
        elif file_type == XLS_FILE and sheet_name != '--':
            df_merge = pd.read_excel(merge_file, sheet_name=sheet_name, skiprows=skip_rows, usecols=lambda x: 'Unnamed' not in x)
            df_merge = df_merge.dropna(how='all')
        else:
            pass
    else:
        pass


    ## dispaly the merge df on top
    # if df_merge is not None:
    #     merge_container.write('External file to be merged')
    #     merge_container.dataframe(df_merge)

    if df_upload is not None and df_merge is not None:
        upload_col, merge_col = st.columns(2)
        upload_headers = list(df_upload)
        upload_headers.insert(0,'--')
        merge_headers = list(df_merge)
        merge_headers.insert(0, '--')
        left_on = upload_col.selectbox('File merged on:', upload_headers, key='merge_left')
        right_on = merge_col.selectbox('Merged file on:', merge_headers, key='merge_right')
        if left_on != '--' and right_on != '--':
            df_upload = df_upload.merge(df_merge, left_on=left_on, right_on=right_on, how='left')



value_mapping_expander = st.sidebar.expander("Value mapping or modification")
with value_mapping_expander:
    num_mapping = st.number_input('Numbers of fields to map:', value=0, key='nf_map')
    if num_mapping:
        for i_mapping in range(num_mapping):
            field = st.selectbox(f'Column {i_mapping+1}:', list(df_upload), key=f'field_{i_mapping}')
            mapping = st.text_area(f"Value mapping {i_mapping+1}:", key=f'mapping_{i_mapping}')
            new_field = st.text_input(f'New field name {i_mapping+1}:')

            if mapping:
                if new_field:
                    df_upload[new_field] = df_upload[field].apply(get_value_from_key, kv_dict=get_rename_dict(mapping))
                else:
                    df_upload[field] = df_upload[field].apply(get_value_from_key, kv_dict=get_rename_dict(mapping))

    MOD_ACTIONS = ['--', 'Prefix', 'Append', 'Multiply', 'Divide', 'Log(e based)', 'Log(10 based)', 'e^', '10^']
    num_mod = st.number_input('Numbers of fields to Modify:', value=0, key='nf_mod')
    if num_mod:
        for i_mod in range(num_mod):
            field = st.selectbox(f'Column {i_mod+1}:', list(df_upload), key=f'field_mod_{i_mod}')
            mod_action = st.selectbox('Modification Action:', MOD_ACTIONS, key=f'mod_action_{i_mod}')
            new_field = st.text_input(f'New field name {i_mod+1}:', key=f'new_mod_field_{i_mod}')

            if mod_action == 'Prefix':
                prefix = st.text_input('Prefix:')
                if new_field:
                    df_upload[new_field] = str(prefix) + df_upload[field].astype(str)
                else:
                    df_upload[field] = str(prefix) + df_upload[field].astype(str)
            elif mod_action == 'Append':
                postfix = st.text_input('Postfix:')
                if new_field:
                    df_upload[new_field] = df_upload[field].astype(str) + str(postfix)
                else:
                    df_upload[field] = df_upload[field].astype(str) + str(postfix)
            elif mod_action == 'Multiply':
                factor = st.number_input("Multiply by:", value=1.00, key='multi')
                if new_field:
                    df_upload[new_field] = df_upload[field].apply(lambda x: float(x) * factor)
                else:
                    df_upload[field] = df_upload[field].apply(lambda x: float(x) * factor)
            elif mod_action == 'Divide':
                factor = st.number_input("Divide by:", value=1.00, key='divide')
                if new_field:
                    df_upload[new_field] = df_upload[field].apply(lambda x: float(x) / factor)
                else:
                    df_upload[field] = df_upload[field].apply(lambda x: float(x) / factor)
            elif mod_action == 'Log(e based)':
                if new_field:
                    df_upload[new_field] = df_upload[field].apply(lambda x: math.log(x))
                else:
                    df_upload[field] = df_upload[field].apply(lambda x: math.log(x))
            elif mod_action == 'Log(10 based)':
                if new_field:
                    df_upload[new_field] = df_upload[field].apply(lambda x: math.log10(x))
                else:
                    df_upload[field] = df_upload[field].apply(lambda x: math.log10(x))
            elif mod_action == 'e^':
                if new_field:
                    df_upload[new_field] = df_upload[field].apply(lambda x: math.exp(x))
                else:
                    df_upload[field] = df_upload[field].apply(lambda x: math.exp(x))
            elif mod_action == '10^':
                if new_field:
                    df_upload[new_field] = df_upload[field].apply(lambda x: math.pow(10, x))
                else:
                    df_upload[field] = df_upload[field].apply(lambda x: math.pow(10, x))

    #
    num_cla = st.number_input('Numbers of fields to classify:', value=0, key='classify')
    if num_cla:
        for i_cla in range(num_cla):
            field = st.selectbox(f'Column {i_cla + 1}:', list(df_upload), key=f'field_cla_{i_cla}')
            class_map = st.text_area(f"Class map {i_cla + 1} (ex. --n2=c1, n2--n3=c2, n3--=c3 )", key=f'cla_{i_cla}')
            new_field = st.text_input(f'New field name {i_cla + 1}:')
            new_field_dub = f'{new_field}_DUB'
            if class_map and new_field:
                df_upload[new_field] = ''
                class_dict = get_rename_dict(class_map)
                for k, v in class_dict.items():
                    ranges = k.split('--')
                    # print('ranges = ', ranges)
                    low = ranges[0].strip()
                    high = ranges[1].strip()
                    if low and high:
                        df_upload[new_field_dub] = df_upload[field].apply(lambda x: v if float(low) < float(x) <= float(high) else '')
                        df_upload[new_field] = df_upload[new_field]+df_upload[new_field_dub]
                    elif high:
                        df_upload[new_field_dub] = df_upload[field].apply(lambda x: v if float(x) <= float(high) else '')
                        df_upload[new_field] = df_upload[new_field] + df_upload[new_field_dub]
                    elif low:
                        df_upload[new_field_dub] = df_upload[field].apply(lambda x: v if float(x) > float(low) else '')
                        df_upload[new_field] = df_upload[new_field] + df_upload[new_field_dub]


                del df_upload[new_field_dub]
    #
    num_nan_fill = st.number_input('Numbers of fields to fill NaN:', value=0, key='fill_nan')
    if num_nan_fill:
        for i_nan_fill in range(num_nan_fill):
            field = st.selectbox(f'Column {i_nan_fill+1}:', list(df_upload), key=f'field_nan_{i_nan_fill}')
            fill = st.text_input(f"Fill NaN/NaT with {i_nan_fill+1}:", key=f'fill_nan_{i_nan_fill}')
            new_field = st.text_input(f'New field name {i_nan_fill+1}:', key=f'fill_new_{i_nan_fill}')

            if fill:
                if new_field:
                    df_upload[new_field] = df_upload[field].apply(lambda x: fill if pd.isna(x) else x)
                else:
                    df_upload[field] = df_upload[field].apply(lambda x: fill if pd.isna(x) else x)


multi_column_expander = st.sidebar.expander("Multi-Columns Operations")
with multi_column_expander:
    headers = list(df_upload)
    headers_select = headers.copy()
    headers_select.insert(0, '--')
    num_merge = st.number_input('Numbers of pairs to merge:', value=0, key='mc_merge')
    if num_merge > 0:
        for group in range(num_merge):
            header_1 = st.selectbox(f'Select header 1 for group {group + 1}', headers_select, key=f'mc_header1{group}')
            header_2 = st.selectbox(f'Select header 2 for group {group + 1}', headers_select, key=f'mc_header2{group}')
            sep = st.text_input('Separated by', value=',', key=f'mc_sep{group}')
            header_new = st.text_input('New header name', value='', key=f'mc_new{group}')
            if header_1 != '--' and header_2 != '--' and header_new:
                df_upload[header_new] = df_upload[header_1].astype(str) + sep + df_upload[header_2].astype(str)

                headers = list(df_upload)
                headers_select = headers.copy()
                headers_select.insert(0, '--')

    num_op = st.number_input('Numbers of pairs to operate on:', value=0, key='mc_op')
    ops = ['+', '-', '*', '/']
    if num_op > 0:
        for group in range(num_op):
            header_1 = st.selectbox(f'Select header 1 for group {group + 1}', headers_select, key=f'mcop_header1{group}')
            header_2 = st.selectbox(f'Select header 2 for group {group + 1}', headers_select, key=f'mcop_header2{group}')
            op = st.selectbox('Operation', ops, key=f'mc_op{group}')

            header_new = st.text_input('New header name', value='', key=f'mcop_new{group}')
            if header_1 != '--' and header_2 != '--' and header_new:
                if op == '+':
                    df_upload[header_new] = df_upload[header_1].astype(float) + df_upload[header_2].astype(float)
                elif op == '-':
                    df_upload[header_new] = df_upload[header_1].astype(float) - df_upload[header_2].astype(float)
                elif op == '*':
                    df_upload[header_new] = df_upload[header_1].astype(float) * df_upload[header_2].astype(float)
                elif op == '/':
                    df_upload[header_new] = df_upload[header_1].astype(float) / df_upload[header_2].astype(float)

                headers = list(df_upload)
                headers_select = headers.copy()
                headers_select.insert(0, '--')

type_expander = st.sidebar.expander("Specify data type")
with type_expander:
    all_headers = list(df_upload.columns.values)
    date_headers = st.multiselect('The following are convert to Datetime type', all_headers, key='date_headers')
    if date_headers and len(date_headers) > 0:
        for col in date_headers:
            df_upload[col] = pd.to_datetime(df_upload[col])
    numeric_headers = st.multiselect('The following are convert to Numerical type', all_headers, key='numeric_headers')
    if numeric_headers and len(numeric_headers) > 0:
        for col in numeric_headers:
            df_upload[col] = pd.to_numeric(df_upload[col])
    string_headers = st.multiselect('The following are convert to String type', all_headers, key='string_headers')
    if string_headers and len(string_headers) > 0:
        df_upload[ string_headers] = df_upload[string_headers].astype(str)


rename_expander = st.sidebar.expander("Rename and Filter")
filter_by = '--'
filter_value = '--'
with rename_expander:
    header_rename = rename_expander.text_area('Rename Headers: (old1=new1,old2=new2, ....)')
    rename_dict = get_rename_dict(header_rename)

    headers = list(df_upload)
    if rename_dict and len(rename_dict) > 0:
        #df_upload = df_upload[headers]
        df_upload = df_upload.rename(columns=rename_dict)
    else:
        df_upload = df_upload[headers]

    st.write('***')
    filter_headers = list(df_upload)
    filter_headers.insert(0, '--')

    num_filter_nn = st.number_input('Number of Not-NA filters', 0, key='filter_nn')
    if num_filter_nn:
        for i_filter_nn in range(num_filter_nn):
            filter_not_none = st.selectbox('Remove Null Field:', filter_headers, key=f'filter_not_none-{i_filter_nn}')
            if filter_not_none != '--':
                df_upload = df_upload[df_upload[filter_not_none].notnull()]

    st.write('***')
    num_filter = st.number_input('Number of value filters', 0, key='nf')
    if num_filter:
        for i_filter in range(num_filter):
            filter_by = st.selectbox('Filter by:', filter_headers, key=f'filter{i_filter}')
            if filter_by != '--':
                data_type_str = str(df_upload[filter_by].dtype)
                #string
                if data_type_str.startswith('object'):
                    part = st.text_input('Partial text:', key=f't_inut_{i_filter}')
                    df_upload['_PART_'] = df_upload[filter_by].apply(
                        lambda t: 'N' if pd.isna(t) else ('Y' if part.lower() in t.lower() else 'N') )
                    df_upload = df_upload[df_upload['_PART_'] == 'Y']
                    del df_upload['_PART_']
                elif data_type_str.startswith('datetime'):
                    filter_date = st.date_input('Input a Date:', key=f'd_inut_{i_filter}')
                    filter_date_str = filter_date.strftime('%Y-%m-%d')
                    #print(filter_date_str)
                    filter_op = st.selectbox("Filter operator:", options=['--', '=', '>', '<', '>=', '<='], key=f'op_date_{i_filter}')
                    if filter_op != '--':
                        if filter_op == '=':
                            df_upload = df_upload.loc[df_upload[filter_by] == filter_date_str]
                        elif filter_op == '>':
                            df_upload = df_upload.loc[df_upload[filter_by] > filter_date_str]
                        elif filter_op == '<':
                            df_upload = df_upload.loc[df_upload[filter_by] < filter_date_str]
                        elif filter_op == '>=':
                            df_upload = df_upload.loc[df_upload[filter_by] >= filter_date_str]
                        elif filter_op == '<=':
                            df_upload = df_upload.loc[df_upload[filter_by] <= filter_date_str]
                elif data_type_str.startswith('int') or data_type_str.startswith('float'):
                    filter_value = st.number_input('Input a Number:', key=f'n_inut_{i_filter}')
                    filter_op = st.selectbox("Filter operator:", options=['--', '=', '>', '<', '>=', '<='], key=f'op_num_{i_filter}')
                    if filter_op != '--':
                        if filter_op == '=':
                            df_upload = df_upload.loc[df_upload[filter_by] == filter_value]
                        elif filter_op == '>':
                            df_upload = df_upload.loc[df_upload[filter_by] > filter_value]
                        elif filter_op == '<':
                            df_upload = df_upload.loc[df_upload[filter_by] < filter_value]
                        elif filter_op == '>=':
                            df_upload = df_upload.loc[df_upload[filter_by] >= filter_value]
                        elif filter_op == '<=':
                            df_upload = df_upload.loc[df_upload[filter_by] <= filter_value]

    st.write('***')
    num_filter_out = st.number_input('Number of filters to exclude from a list.', 0, key='exclude')
    if num_filter_out:
        for i_filter_out in range(num_filter_out):
            filter_out_by = st.selectbox('The following are excluded from:', filter_headers,
                                         key=f'filter_out_by_{i_filter_out}')
            if filter_out_by != '--':
                filter_out_string = st.text_area(r'List to exclude (separate by , or newline:',
                                               key=f'filter_out_{i_filter_out}')
                filter_out_list = get_list(filter_out_string)
                data_type_str = str(df_upload.dtypes[filter_out_by])
                if data_type_str.startswith('object') or data_type_str.startswith('int'):
                    if filter_out_list:
                        df_upload['_PART_'] = df_upload[filter_out_by].apply(lambda t: 'Y' if pd.isna(t) else ('Y' if str(t) not in filter_out_list else 'N'))
                        df_upload = df_upload[df_upload['_PART_'] == 'Y']
                        del df_upload['_PART_']
                else:
                    st.info('Not support for list exclusion.')

    st.write('***')
    num_filter_in = st.number_input('Number of filters to include only from a list.', 0, key='include')
    if num_filter_in:
        for i_filter_in in range(num_filter_in):
            filter_in_by = st.selectbox('The following are only included from:', filter_headers, key=f'filter_in_by_{i_filter_in}')
            if filter_in_by != '--':
                filter_in_string = st.text_area(r'List to include (separate by , or newline:', key=f'filter_in_{i_filter_in}')
                filter_in_list = get_list(filter_in_string)
                data_type_str = str(df_upload.dtypes[filter_in_by])
                if data_type_str.startswith('object') or data_type_str.startswith('int'):
                    if filter_in_list:
                        df_upload['_PART_'] = df_upload[filter_in_by].apply(lambda t: 'N' if pd.isna(t) else ('Y' if str(t) in filter_in_list else 'N'))
                        df_upload = df_upload[df_upload['_PART_'] == 'Y']
                        del df_upload['_PART_']
                else:
                    st.info('Not support for list only inclusion.')


headers = list(df_upload)


add_cfield_expander = st.sidebar.expander("Add constant fields")
with add_cfield_expander:
    column_value = st.text_area('Column names with values: (col1=val1,col2=val2, ....)', key='col_val')
    column_value_dict = get_rename_dict(column_value)

    if column_value_dict and len(column_value_dict) > 0:
        for key in column_value_dict:
            df_upload[key] = column_value_dict[key]


headers = list(df_upload)

fields_expander = st.sidebar.expander("Select fields and sorting")
with fields_expander:
    all_headers = list(df_upload.columns.values)
    headers = fields_expander.multiselect('Select headers', all_headers, all_headers)
    df_upload = df_upload[headers]

    sort_result = st.radio('Sort Results?', ['No', 'Ascending', 'Decending'])
    if sort_result == 'Ascending' or sort_result == 'Decending':
        sort_headers = ['--'] + headers.copy()
        # sort_headers.remove(STRUCTURE)
        sort_by = fields_expander.selectbox('Sort by:', sort_headers)
        if sort_by != '--':
            ASCENDING = sort_result == 'Ascending'
            df_upload.sort_values(by=sort_by, ascending=ASCENDING, inplace=True)


headers = list(df_upload)

df_upload = df_upload[headers]
if STRUCTURE in df_upload:
    df_upload.drop(STRUCTURE, axis=1, inplace=True)

st.write(f'Data size = {len(df_upload)}')
event = st.dataframe(
        df_upload,
        key="df",
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
    )

if SMILES in df_upload:
    selected_rows = event.selection['rows']
    for row in selected_rows:
        smi = df_upload.at[row, SMILES]
        mol = Chem.MolFromSmiles(smi)
        details_expander.write( moltosvg(mol), unsafe_allow_html=True) 
    


if df_upload is not None:
    date_str = str(date.today())

    with excel_download_container.popover('Download Xlsx file'):
        sheet_name = st.text_input("Sheet name:")
        excel_download_file_name = f'{upload_file_name}_{date_str}.xlsx'
        st.download_button('Download Xlsx file', get_df_xlsx(df_upload, sheet_name), excel_download_file_name)
    # smaller_size = download_container.number_input('Download Xlsx in smaller size of:', 0, step=500, key='smaller')
    # if smaller_size:
    #     df_chunk_list = split_dataframe(df_upload, smaller_size)
    #     for idx, df_trunk in enumerate(df_chunk_list):
    #         download_file_name = f'{upload_file_name}_{date_str}_{idx+1}.xlsx'
    #         download_container.download_button(f'Download [{idx+1}]', get_df_xlsx(df_trunk), download_file_name)


    with csv_download_container.popover('Download csv file'):
        delimiter = st.selectbox('Delimiter', options=[',', 'tab', ';'])
        if delimiter == 'tab':
            delimiter = '\t'
        extension = st.text_input('extension', value='csv')
        csv_download_file_name = f'{upload_file_name}_{date_str}.{extension}'
        st.download_button('Download CSV file', get_df_csv(df_upload, sep=delimiter), csv_download_file_name)
    # csv_smaller_size = csv_download_container.number_input('Download CSV in smaller size of:', 0, step=500, key='csv_smaller')
    # if csv_smaller_size:
    #     df_chunk_list = split_dataframe(df_upload, csv_smaller_size)
    #     for idx, df_trunk in enumerate(df_chunk_list):
    #         csv_download_file_name = f'{upload_file_name}_{date_str}_{idx+1}.csv'
    #         csv_download_container.download_button(f'Download [{idx+1}]', get_df_csv(df_trunk, sep=delimiter), csv_download_file_name)

    sdf_download_file_name = f'{upload_file_name}_{date_str}.sdf'
    if SMILES in df_upload:
        PandasTools.AddMoleculeColumnToFrame(df_upload, SMILES, STRUCTURE, includeFingerprints=False)
        sdf_download_container.download_button('Download Sdf file', get_df_sdf(df_upload, STRUCTURE), sdf_download_file_name)
    
    # df_headers = pd.DataFrame({'Header': list(df_upload)})
    # headers_container.dataframe(df_headers)