import sys

def set_column_format(table, name, format):
    if name in table.colnames:
        table[name].info.format = format

def default_coord_format(c):
    return c.to_string(style='dms', format='unicode', precision=1)

def default_mag_err_format(m):
    return f"{m.value[0]:.2f} ± {m.value[1]:.3f}"

def default_table_format(table):
    set_column_format(table, 'flux', '.0f')
    set_column_format(table, 'snr', '.1f')
    set_column_format(table, 'peak', '.1%')
    set_column_format(table, 'radec2000', default_coord_format)
    set_column_format(table, 'sky_centroid', default_coord_format)
    set_column_format(table, 'M', default_mag_err_format)
