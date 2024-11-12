import re
import gzip

import hail as hl

from typing import Union

# ensure that PYTHONPATH includes the wes-qc directory
from utils.config import path_spark


def compare_structs(struct1, struct2):
    """
    Compare each field in two struct entries.
    
    Returns:
    bool: True if entries are equal, False otherwise.
    """
    fields = list(struct1.dtype.fields)

    # if not fields other than key
    if len(fields) == 0:
        return hl.bool(True)
    
    # Check for NA fields
    def check_field_equal(field):
        # both fields not defined
        field_not_defined = (~hl.is_defined(struct1[field]) & ~hl.is_defined(struct2[field]))

        # both fields defined
        field_defined = hl.is_defined(struct1[field]) & hl.is_defined(struct2[field])

        # different checks for differenet field dtypes
        if struct1[field].dtype == hl.dtype('float'):
            # both fields not a number
            field_nan = hl.is_nan(struct1[field]) & hl.is_nan(struct2[field])
            expr = (field_defined & hl.approx_equal(struct1[field], struct2[field])) | field_not_defined | field_nan
        elif struct1[field].dtype == hl.dtype('array<float>'):
            # element_wise_expr = hl.map(lambda a, b: (hl.is_defined(a) & hl.is_defined(b) & hl.approx_equal(a, b)) | (~hl.is_defined(a) & ~hl.is_defined(b)), struct1[field], struct2[field])

            element_wise_expr = hl.map(lambda a, b: (hl.is_defined(a) &
                                                     hl.is_defined(b) &
                                                     hl.approx_equal(a, b)) | 
                                                     (~hl.is_defined(a) & 
                                                      ~hl.is_defined(b)) |
                                                     (hl.is_nan(a) &
                                                      hl.is_nan(b)),
                                       struct1[field], struct2[field])
            expr = hl.all(element_wise_expr)
        else:
            # TODO: implement comparison for other dtypes
            expr = (field_defined & (struct1[field] == struct2[field])) | field_not_defined
        
        return expr

    comparisons = [check_field_equal(field) for field in fields]
    return hl.all(lambda x: x, comparisons)

def compare_entries_to_other_mt(this_entry, row_key, col_key, other_mt):
    """
    Compare entries of two MatrixTables.
    
    Returns:
    bool: True if entries are equal, False otherwise.
    """
    other_entry = other_mt[row_key, col_key]
    return ~compare_structs(this_entry, other_entry)

def compare_plinks(bed1_path: str, bim1_path: str, fam1_path: str, 
                   bed2_path: str, bim2_path: str, fam2_path: str) -> bool:
    """
    Compare the contents of two Hail plink files.

    Returns:
    bool: True if the plinks are equal, False otherwise.
    """
    # Load the plinks
    mt1 = hl.import_plink(bed=path_spark(bed1_path), bim=path_spark(bim1_path), fam=path_spark(fam1_path), reference_genome='GRCh38')
    mt2 = hl.import_plink(bed=path_spark(bed2_path), bim=path_spark(bim2_path), fam=path_spark(fam2_path), reference_genome='GRCh38')

    # Ensure the schemas are the same
    if mt1.row.dtype != mt2.row.dtype or mt1.col.dtype != mt2.col.dtype or mt1.entry.dtype != mt2.entry.dtype:
        print('MatrixTable schemas do not match')
        return False

    # Align row and column keys
    mt2 = mt2.key_rows_by(*mt1.row_key)  # Reorder row keys to match mt1
    mt2 = mt2.key_cols_by(*mt1.col_key)  # Reorder column keys to match mt1
    
    filtered_mt = mt1.filter_entries(compare_entries_to_other_mt(mt1.entry, mt1.row_key, mt1.col_key, mt2))

    # Check if there are any differing entries
    num_differences = filtered_mt.entries().count()
    if num_differences == 0:
        print('MatrixTables are equal')
        return True
    else:
        print(f'MatrixTables are not equal: {num_differences} differing entries found')
        return False

def compare_matrixtables(mt1_path: Union[str, hl.MatrixTable], mt2_path: Union[str, hl.MatrixTable]) -> bool:
    """
    Compare the contents of two Hail MatrixTables.

    Parameters:
    mt1_path (str): Path to the first MatrixTable.
    mt2_path (str): Path to the second MatrixTable.

    Returns:
    bool: True if the MatrixTables are equal, False otherwise.
    """
    # Load the MatrixTables
    if isinstance(mt1_path, str):
        mt1 = hl.read_matrix_table(path_spark(mt1_path))
    elif isinstance(mt1_path, hl.MatrixTable):
        mt1 = mt1_path
    else:
        raise TypeError(f"Expected 'str' or 'hl.Table', but got {type(mt1_path)}")
    
    if isinstance(mt2_path, str):
        mt2 = hl.read_matrix_table(path_spark(mt2_path))
    elif isinstance(mt2_path, hl.MatrixTable):
        mt2 = mt2_path
    else:
        raise TypeError(f"Expected 'str' or 'hl.Table', but got {type(mt2_path)}")

    # Ensure the schemas are the same
    if mt1.row.dtype != mt2.row.dtype or mt1.col.dtype != mt2.col.dtype or mt1.entry.dtype != mt2.entry.dtype:
        print('MatrixTable schemas do not match')
        return False

    # Align row and column keys
    mt2 = mt2.key_rows_by(*mt1.row_key)  # Reorder row keys to match mt1
    mt2 = mt2.key_cols_by(*mt1.col_key)  # Reorder column keys to match mt1
    
    filtered_mt = mt1.filter_entries(compare_entries_to_other_mt(mt1.entry, mt1.row_key, mt1.col_key, mt2))

    # Check if there are any differing entries
    num_differences = filtered_mt.entries().count()
    if num_differences == 0:
        print('MatrixTables are equal')
        return True
    else:
        print(f'MatrixTables {mt1_path}\n{mt2_path}\nare not equal: {num_differences} differing rows found')
        return False

def compare_entries_to_other_ht(this_row_value, row_key, other_ht):
    """
    Compare rows of two Tables.
    
    Returns:
    bool: True if entries are equal, False otherwise.
    """
    other_row_value = other_ht[row_key]
    return ~compare_structs(this_row_value, other_row_value)

def compare_tables(ht1_path: Union[str, hl.Table], ht2_path: Union[str, hl.Table]) -> bool:
    """
    Compare the contents of two Hail Tables.

    Parameters:
    ht1_path (str): Path to the first Table.
    ht2_path (str): Path to the second Table.

    Returns:
    bool: True if the Tables are equal, False otherwise.
    """
    # Load the Tables
    if isinstance(ht1_path, str):
        ht1 = hl.read_table(path_spark(ht1_path))
    elif isinstance(ht1_path, hl.Table):
        ht1 = ht1_path
    else:
        raise TypeError(f"Expected 'str' or 'hl.Table', but got {type(ht1_path)}")
    
    if isinstance(ht2_path, str):
        ht2 = hl.read_table(path_spark(ht2_path))
    elif isinstance(ht2_path, hl.Table):
        ht2 = ht2_path
    else:
        raise TypeError(f"Expected 'str' or 'hl.Table', but got {type(ht2_path)}")

    # Ensure the schemas are the same
    if ht1.row.dtype != ht2.row.dtype:
        print("Table schemas do not match")
        return False

    # Align row keys
    ht2 = ht2.key_by(*ht1.key)  # Reorder row keys to match ht1

    rows_differ = ht1.filter(compare_entries_to_other_ht(ht1.row_value, ht1.key, ht2))

    # Check if there are any differing rows
    num_differences = rows_differ.count()

    if num_differences == 0:
        print(f'Tables {ht1_path}\n{ht2_path}\nare equal')
        return True
    else:
        print(f'Tables {ht1_path}\n{ht2_path}\nare not equal: {num_differences} differing rows found')
        return False

def compare_txts(path_1: str, path_2: str, replace_strings: list[list[str, str]] = None) -> bool:
    """
    Compare contents of two files based on their raw string content.

    Returns:
    bool: True if files' contents are equal, False otherwise.
    """
    with open(path_1, 'r') as f_1:
        contents_1 = f_1.read()

    with open(path_2, 'r') as f_2:
        contents_2 = f_2.read()

    if replace_strings:
        for pattern, substitute in replace_strings:
            contents_1 = re.sub(pattern, substitute, contents_1)
            contents_2 = re.sub(pattern, substitute, contents_2)

    if contents_1 == contents_2:
        print(f'Files {path_1}\n{path_2}\n are equal')
    else:
        print(f'Files {path_1}\n{path_2}\n are not equal')

    return contents_1 == contents_2

def compare_bgzed_txts(path_1: str, path_2: str, replace_strings: list[list[str, str]] = None) -> bool:
    """
    Compare contents of two bgzed txt files. 
    
    Parameters:
    replace_strings: list[list[str, str]] - regex and str to replace

    Returns:
    bool: True if files' contents are equal, False otherwise.
    """
    with gzip.open(path_1, 'r') as f_1:
        contents_1 = f_1.read()

    with gzip.open(path_2, 'r') as f_2:
        contents_2 = f_2.read()
    
    if replace_strings:
        for pattern, substitute in replace_strings:
            contents_1 = re.sub(pattern, substitute, contents_1)
            contents_2 = re.sub(pattern, substitute, contents_2)

    if contents_1 == contents_2:
        print(f'Files {path_1}\n{path_2}\n are equal')
    else:
        print(f'Files {path_1}\n{path_2}\n are not equal')

    return contents_1 == contents_2
