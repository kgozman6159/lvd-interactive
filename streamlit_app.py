# ----------import necessary packages---------- #
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import altair as alt

from astropy.io import fits
import astropy.table as table
from astropy.io import ascii

from astropy import units as u

import astropy.coordinates as coord

import local_volume_database # idk what this does yet

#st.set_page_config(layout="centered")
# ---------------------load data---------------------- #
@st.cache_data
def load_data():
    # loads versions from latest github release
    dwarf_all = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_all.csv')
    dsph_mw = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_mw.csv')
    dsph_m31 = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_m31.csv')
    dsph_lf = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_local_field.csv')
    dsph_lf_distant = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_local_field_distant.csv')
    gc_ambiguous = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/gc_ambiguous.csv')
    gc_mw_new = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/gc_mw_new.csv')
    gc_harris = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/gc_harris.csv')
    gc_dwarf_hosted = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/gc_dwarf_hosted.csv')
    gc_other = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/gc_other.csv')
    candidate = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/candidate.csv')
    misc_host = pd.read_csv('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/misc_host.csv')

    # ----don't know if I need thes below since dsph_mw already has columns for the wolf dynamical mass and HI mass UL which are very similar ----#

    #dsph_mw['mass_dynamical'] = dsph_mw['vlos_sigma']**2*930*dsph_mw['rhalf']*np.sqrt(1.-dsph_mw['ellipticity'])* dsph_mw['distance']*np.pi/180./60.*1000.
    #dsph_mw['mass_dynamical_ul'] = dsph_mw['vlos_sigma_ul']**2*930*dsph_mw['rhalf']*np.sqrt(1.-dsph_mw['ellipticity'])* dsph_mw['distance']*np.pi/180./60.*1000.

    #dsph_mw['mass_HI_ul'] = (235600*dsph_mw['flux_HI_ul']*(dsph_mw['distance']/1000.)**2) # this is just the 10**x version of the mass_HI_ul column in dwarf_all
    #comb['mass_HI_ul'] = np.log10(235600*comb['flux_HI_ul']*(comb['distance']/1000.)**2)

    return dwarf_all, dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate, misc_host

dwarf_all, dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate, misc_host = load_data()


# ---------------------dictionary of column labels and descriptions---------------------- #

table_descriptions = pd.read_csv('table_descriptions.csv', index_col='property', keep_default_na=False)
#print(table_descriptions.iloc[0].unit)
table_descriptions = table_descriptions.T.to_dict(index='property')
print(table_descriptions.keys())
#st.write(table_descriptions)
# st.write(table_descriptions['rhalf_physical']['desc'])
valid_plot_cols = ['key', 'ra', 'dec', 'name', 'host', 'confirmed_real', 
                   'confirmed_dwarf', 'rhalf',  'position_angle', 'ellipticity', 
                   'ellipticity_ul', 'distance_modulus', 'apparent_magnitude_v', 
                   'vlos_systemic', 'vlos_sigma', 'vlos_sigma_ul', 'pmra',  'pmdec', 'metallicity_spectroscopic', 
                'metallicity_spectroscopic_sigma',  'metallicity_spectroscopic_sigma_ul', 'rcore', 'rking', 'rad_sersic', 
                'n_sersic', 'age', 'metallicity_isochrone', 'flux_HI', 'flux_HI_ul', 'metallicity_photometric', 
                'metallicity_photometric_sigma', 'metallicity_photometric_sigma_ul', 'M_V', 'mass_stellar', 'distance',
                'll', 'bb', 'sg_xx', 'sg_yy', 'sg_zz', 'distance_gc', 'distance_m31', 'distance_lg', 'distance_host', 
                'mass_HI', 'mass_HI_ul', 'metallicity', 'metallicity_type', 'velocity_gsr', 'velocity_lg', 
                'mass_dynamical_wolf', 'mass_dynamical_wolf_ul', 'rhalf_physical', 'rhalf_sph_physical', 'surface_brightness_rhalf']

# ---------------------misc. functions---------------------- #
# ------ M_V <-> L_V ------ #
@st.cache_data
def lum(x):
    m_x_sun=4.83
    return -.4*(x - m_x_sun) + np.log10(2.)

@st.cache_data
def lum_inverse(x):
    m_x_sun=4.83
    return m_x_sun - (x )/0.4 + np.log10(2.)

#type_x, type_y = 'linear', 'linear'
#reverse_y, reverse_x = False, False


def get_axis_specs(axis, key):
    type_axis = 'linear'
    if table_descriptions[axis]['dtype'] in ['float64']:
        if not (dwarf_all[axis] < 0).any():
            type_axis = st.sidebar.pills(axis + ' scale', ['linear', 'log'], default='linear', key=key)

    if axis in ['apparent_magnitude_v', "M_V"]:
        reverse_axis = True
    else:
        reverse_axis = False
    if table_descriptions[axis]['unit'] != "N/A":
        axis_label = table_descriptions[axis]['label'] + ' (' + table_descriptions[axis]['unit'] + ')'
    else:
        axis_label = table_descriptions[axis]['label']

    return type_axis, reverse_axis, axis_label



plot_xaxis = st.sidebar.selectbox('x-axis', valid_plot_cols)
type_x, reverse_x, xlabel = get_axis_specs(plot_xaxis, 'xaxis')

# if table_descriptions[plot_xaxis]['dtype'] in ['float64']:
#     if not (dwarf_all[plot_xaxis] < 0).any():
#         type_x = st.sidebar.pills('x-axis scale', ['linear', 'log'], default='linear', key='xscale')


plot_yaxis = st.sidebar.selectbox('y-axis', valid_plot_cols)
type_y, reverse_y, ylabel = get_axis_specs(plot_yaxis, 'yaxis')

# if table_descriptions[plot_yaxis]['dtype'] in ['float64']:
#     if not (dwarf_all[plot_yaxis] < 0).any():
#         type_y = st.sidebar.pills('y-axis scale', ['linear', 'log'], default='linear', key='yscale')

# if plot_xaxis in ['apparent_magnitude_v', "M_V"]:
#     reverse_x = True
# if plot_yaxis in ['apparent_magnitude_v', "M_V"]:
#     reverse_y = True

# if table_descriptions[plot_xaxis]['unit'] != "N/A":
#     xlabel = table_descriptions[plot_xaxis]['label'] + ' (' + table_descriptions[plot_xaxis]['unit'] + ')'
# else:
#     xlabel = table_descriptions[plot_xaxis]['label']

plot = st.altair_chart(alt.Chart(dsph_mw).mark_circle().encode(
     x=alt.X(plot_xaxis, scale=alt.Scale(type=type_x, reverse=reverse_x), title=xlabel), 
     y=alt.Y(plot_yaxis, scale=alt.Scale(type=type_y, reverse=reverse_y), title=table_descriptions[plot_yaxis]['label']), 
     color='host').interactive(), use_container_width=True)


plot = st.empty()

#st.toast('Welcome to the Local Volume Database!')





st.text('all dwarfs')
st.dataframe(dwarf_all, use_container_width=True) # use_container_width doesn't work??
st.text('dSphs in MW')
st.dataframe(dsph_mw, use_container_width=False)
#st.dataframe(misc_host)






df2 = pd.DataFrame(
     np.random.randn(200, 3),
     columns=['a', 'b', 'c'])
c = alt.Chart(df2).mark_circle().encode(
     x='a', y='b', size='c', color='c', tooltip=['a', 'b', 'c']).interactive()
#st.write(c)