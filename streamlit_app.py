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

    dsph_mw['mass_dynamical'] = dsph_mw['vlos_sigma']**2*930*dsph_mw['rhalf']*np.sqrt(1.-dsph_mw['ellipticity'])* dsph_mw['distance']*np.pi/180./60.*1000.
    dsph_mw['mass_dynamical_ul'] = dsph_mw['vlos_sigma_ul']**2*930*dsph_mw['rhalf']*np.sqrt(1.-dsph_mw['ellipticity'])* dsph_mw['distance']*np.pi/180./60.*1000.

    dsph_mw['mass_HI_ul'] = (235600*dsph_mw['flux_HI_ul']*(dsph_mw['distance']/1000.)**2)
    #comb['mass_HI_ul'] = np.log10(235600*comb['flux_HI_ul']*(comb['distance']/1000.)**2)

    return dwarf_all, dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate, misc_host

dwarf_all, dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate, misc_host = load_data()

## M_V -> L_V
@st.cache_data
def lum(x):
    m_x_sun=4.83
    return -.4*(x - m_x_sun) + np.log10(2.)
@st.cache_data
def lum_inverse(x):
    m_x_sun=4.83
    return m_x_sun - (x )/0.4 + np.log10(2.)

st.dataframe(dwarf_all, use_container_width=True)
st.dataframe(dsph_mw, use_container_width=False)
#st.dataframe(misc_host)
st.sidebar.write('This is a sidebar')
if st.sidebar.button('Click me'):
    st.sidebar.write('Button clicked!')
if st.button('Say hello'):
     st.write('Why hello there')
else:
     st.write('Goodbye')



st.header('st.write')

# Example 1

st.write('Hello, *World!* :sunglasses:')

# Example 2

st.write(1234)

# Example 3

df = pd.DataFrame({
     'first column': [1, 2, 3, 4],
     'second column': [10, 20, 30, 40]
     })
st.write(df)

# Example 4

st.write('Below is a DataFrame:', df, 'Above is a dataframe.')

# Example 5

df2 = pd.DataFrame(
     np.random.randn(200, 3),
     columns=['a', 'b', 'c'])
c = alt.Chart(df2).mark_circle().encode(
     x='a', y='b', size='c', color='c', tooltip=['a', 'b', 'c']).interactive()
st.write(c)