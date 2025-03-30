import sys
from collections import Counter

def set_column_format(table, name, format):
    if name in table.colnames:
        table[name].info.format = format

def default_coord_format(c):
    return c.to_string(style='dms', format='unicode', precision=1)

def default_mag_err_format(m):
    return f"{m.value[0]:.2f} ± {m.value[1]:.3f}"

def default_val_err_format(m):
    return f"{m.value[0]:.3f} ± {m.value[1]:.3f}"

def default_xfm_val_err_format(m):
    return f"Ta: {m[0][0]:.3f} ± {m[0][1]:.3f}; Tb: {m[1][0]:.3f} ± {m[1][1]:.3f}; Tab: {m[2][0]:.3f} ± {m[2][1]:.3f}"

def default_cxfm_val_err_format(m):
    return f"Ta: {m[0][0]:.3f} ± {m[0][1]:.3f}; Za: {m[1][0]:.3f} ± {m[1][1]:.3f}; Tb: {m[2][0]:.3f} ± {m[2][1]:.3f}; Zb: {m[3][0]:.3f} ± {m[3][1]:.3f}"

def default_table_format(table):
    set_column_format(table, 'flux', '.0f')
    set_column_format(table, 'snr', '.1f')
    set_column_format(table, 'peak', '.1%')
    set_column_format(table, 'radec2000', default_coord_format)
    set_column_format(table, 'sky_centroid', default_coord_format)
    for name in table.colnames:
        column = table.columns[name]
        if hasattr(column, 'dtype'):
            if Counter(column.dtype.names) == Counter(['mag', 'err']):
                set_column_format(table, name, default_mag_err_format)
            elif Counter(column.dtype.names) == Counter(['val', 'err']):
                set_column_format(table, name, default_val_err_format)
            elif Counter(column.dtype.names) == Counter(['Ta', 'Tb', 'Tab']):
                set_column_format(table, name, default_xfm_val_err_format)
            elif Counter(column.dtype.names) == Counter(['Ta', 'Za', 'Tb', 'Zb']):
                set_column_format(table, name, default_cxfm_val_err_format)
    return table
