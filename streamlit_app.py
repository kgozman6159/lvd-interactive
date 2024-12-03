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
from collections import OrderedDict
st.set_page_config(layout="wide")

source = pd.DataFrame({
    "yield_error": [7.5522, 6.9775, 3.9167, 11.9732],
    "yield_center": [32.4, 30.96667, 33.966665, 30.45],
    "variety": ["Glabron", "Manchuria", "No. 457", "No. 462"],
})

bar = alt.Chart(source).mark_errorbar().encode(
    x=alt.X("yield_center:Q").scale(zero=False).title("yield"),
    xError=("yield_error:Q"),
    y=alt.Y("variety:N"),
).transform_filter(
        (alt.datum["yield_error"] != 0))

point = alt.Chart(source).mark_point(
    filled=True,
    color="green",
    opacity=1
).encode(
    alt.X("yield_center:Q"),
    alt.Y("variety:N"),
)

st.altair_chart(alt.layer(bar, point))



table_names = ['dsph_mw', 'dsph_m31', 'dsph_lf', 'dsph_lf_distant', 'gc_ambiguous', 'gc_mw_new', 'gc_harris', 'gc_dwarf_hosted', 'gc_other', 'candidate']
table_names_pretty = ['MW Dwarfs', "M31 Dwarfs", 'Local Field Dwarfs', 'Distant Local Field Dwarfs', 'Ambiguous GCs', 'New MW GCs', 'Harris GCs', 'Dwarf Hosted GCs', 'Other GCs', 'Candidates']

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
    # Combine all tables except dwarf_all into one big dataframe
    tables = [dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate]
    combined_df = pd.concat([table.assign(source=name, source_pretty=name_pretty) for table, name, name_pretty in zip(tables, table_names, table_names_pretty)], ignore_index=True)

    for key in combined_df.keys():
        if (key+'_em' in combined_df.keys()):
            combined_df[key+"_low"] = combined_df[key]-combined_df[key+'_em']
            combined_df[key+"_low"] = combined_df[key+"_low"].fillna(0)
            # print(np.min(combined_df[key]))
            
            # in rare cases, the lower error causes a value that should be positive (like sersic index or vlos,systemic) to be negative (ie if author reports symmetric errorbars)
            # in these cases, need to set these lower error to 0 so that log plots don't break when errorbars are displayed
            if np.min(combined_df[key])>0 and (combined_df[key+"_low"]<0).any():
                combined_df.loc[combined_df[key+"_low"] < 0, key+"_low"] = 0
        if (key+'_ep' in combined_df.keys()):
            combined_df[key+"_high"] = combined_df[key]+combined_df[key+'_ep']
            combined_df[key+"_high"] = combined_df[key+"_high"].fillna(0)
        if (key+'_ul' in combined_df.keys()):
            # print(32*np.nanstd(dwarf_all[key]))
            combined_df[key+"_upper"] = np.ones(len(combined_df[key+'_ul']))*1000*np.nanstd(combined_df[key])
            combined_df[key+"_upper"] = combined_df[key+"_upper"].fillna(0)

        
        # for key in dwarf_all.keys():
        #     if key.endswith('_em') or key.endswith('_ep'):
        #         dwarf_all[key].fillna(0, inplace=True)

    return dwarf_all, dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate, misc_host, combined_df

dwarf_all, dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate, misc_host, master_df = load_data()
print(len(dwarf_all),len(master_df), len(dsph_mw)+len(dsph_m31)+len(dsph_lf)+len(dsph_lf_distant)+len(gc_ambiguous)+len(gc_mw_new)+len(gc_harris)+len(gc_dwarf_hosted)+len(gc_other)+len(candidate),len(misc_host))
#st.dataframe(dwarf_all, use_container_width=True) # use_container_width doesn't work??
#st.dataframe(master_df, use_container_width=True)
# ---------------------dictionary of column labels and descriptions---------------------- #
#print(dwarf_all['M_V_high'])
tab_desc = pd.read_csv('table_descriptions.csv', index_col='property', keep_default_na=False)
#print(tab_desc.iloc[0].unit)
tab_desc = tab_desc.T.to_dict(index='property')
# print(tab_desc.keys())
#st.write(tab_desc)
# st.write(tab_desc['rhalf_physical']['desc'])
valid_plot_cols = ['key', 'ra', 'dec', 'name', 'host', 'confirmed_real', 
                   'confirmed_dwarf', 'rhalf',  'position_angle', 'ellipticity', 
                   'distance_modulus', 'apparent_magnitude_v', 
                   'vlos_systemic', 'vlos_sigma', 'pmra',  'pmdec', 'metallicity_spectroscopic', 
                'metallicity_spectroscopic_sigma', 'rcore', 'rking', 'rad_sersic', 
                'n_sersic', 'age', 'metallicity_isochrone', 'flux_HI', 'metallicity_photometric', 
                'metallicity_photometric_sigma',  'M_V', 'mass_stellar', 'distance',
                'll', 'bb', 'sg_xx', 'sg_yy', 'sg_zz', 'distance_gc', 'distance_m31', 'distance_lg', 'distance_host', 
                'mass_HI', 'metallicity', 'metallicity_type', 'velocity_gsr', 'velocity_lg', 
                'mass_dynamical_wolf', 'rhalf_physical', 'rhalf_sph_physical', 'surface_brightness_rhalf']

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




def get_axis_specs(axis, key):
    type_axis = 'linear'
    reverse_axis = False
    label = tab_desc[axis]['label']
    if tab_desc[axis]['dtype'] in ['float64'] and axis not in ['ra', 'dec']:
        if not (master_df[axis] <= 0).any():
            type_axis = st.segmented_control(label + ' scale', ['linear', 'log'], default='linear', key=key)
    if tab_desc[axis]['dtype'] in ['float64']:
        channel = 'quantitative'
    elif tab_desc[axis]['dtype'] in ['int64']:
        channel = 'ordinal'
    else:
        channel = 'nominal'
    if axis in ['apparent_magnitude_v', "M_V"]:
        reverse_axis = True
    if tab_desc[axis]['unit'] != "N/A":
        axis_label = label + ' (' + tab_desc[axis]['unit'] + ')'
    else:
        axis_label = label
    if type_axis is None:
        type_axis = 'linear'

    if axis+"_high" in master_df.keys():
        show_error = st.checkbox(label + ' error bars?', key=key+"err")
    else:
        show_error = False

    return type_axis, reverse_axis, axis_label, channel, show_error

# ---------------------sidebar---------------------- #

# with st.sidebar:
#     with st.container(border=True, key='xcont') as xcont:
#         left, right = st.columns([7,1], vertical_alignment="top")
#         plot_xaxis = left.selectbox('x-axis', valid_plot_cols, index=1, format_func=lambda x: tab_desc[x]['label'], label_visibility='visible')
#         right.caption("---", help=tab_desc[plot_xaxis]['desc'])
#         type_x, reverse_x, xlabel, channel_x, show_xerr = get_axis_specs(plot_xaxis, 'xaxis')
#     with st.container(border=True, key='ycont') as ycont:
#         left, right = st.columns([7,1], vertical_alignment="top")
#         plot_yaxis = left.selectbox('y-axis', valid_plot_cols, index=2, format_func=lambda x: tab_desc[x]['label'])
#         right.caption("---", help=tab_desc[plot_yaxis]['desc'])
#         type_y, reverse_y, ylabel, channel_y, show_yerr = get_axis_specs(plot_yaxis, 'yaxis')

if "my_key1" not in st.session_state:
    st.session_state.my_key1 = "ra"
if "my_key2" not in st.session_state:
    st.session_state.my_key2 = "dec"

st.session_state.my_key1=st.session_state.my_key1
st.session_state.my_key2=st.session_state.my_key2
with st.sidebar:
    with st.container(border=True, key='xcont') as xcont:
        plot_xaxis = st.selectbox('x-axis', valid_plot_cols, index=1, key="my_key1", format_func=lambda x: tab_desc[x]['label'], label_visibility='visible', help=f"{tab_desc[st.session_state.my_key1]['desc']}")
        #right.caption("---", help=tab_desc[plot_xaxis]['desc'])
        type_x, reverse_x, xlabel, channel_x, show_xerr = get_axis_specs(plot_xaxis, 'xaxis')
    with st.container(border=True, key='ycont') as ycont:
        plot_yaxis = st.selectbox('y-axis', valid_plot_cols, index=2, key="my_key2", format_func=lambda x: tab_desc[x]['label'], help=f"{tab_desc[st.session_state.my_key2]['desc']}")
        #right.caption("---", help=tab_desc[plot_yaxis]['desc'])
        type_y, reverse_y, ylabel, channel_y, show_yerr = get_axis_specs(plot_yaxis, 'yaxis')
# st.selectbox(
#     "Make a selection",
#     ["Default", "Apple", "Banana", "Carrot"],
#     key="my_key",
#     help=f"You've selected {st.session_state.my_key}"
# )


# ---------------------filtering---------------------- #
#filter by source
source = st.sidebar.multiselect('Source', table_names_pretty, default=table_names_pretty)
if source:
    master_df = master_df[master_df['source_pretty'].isin(source)]
# # filter by confirmed dwarf
# confirmed_dwarf = st.sidebar.multiselect('Confirmed Dwarf', dwarf_all['confirmed_dwarf'].unique())
# if confirmed_dwarf:
#     dwarf_all = dwarf_all[dwarf_all['confirmed_dwarf'].isin(confirmed_dwarf)]
# # filter by confirmed real
# confirmed_real = st.sidebar.multiselect('Confirmed Real', dwarf_all['confirmed_real'].unique())
# if confirmed_real:
#     dwarf_all = dwarf_all[dwarf_all['confirmed_real'].isin(confirmed_real)]
# # filter by distance
# distance = st.sidebar.slider('Distance (kpc)', dwarf_all['distance'].min(), dwarf_all['distance'].max(), (dwarf_all['distance'].min(), dwarf_all['distance'].max()))
# dwarf_all = dwarf_all[(dwarf_all['distance'] >= distance[0]) & (dwarf_all['distance'] <= distance[1])]
# # filter by apparent magnitude
# apparent_magnitude_v = st.sidebar.slider('Apparent Magnitude V', dwarf_all['apparent_magnitude_v'].min(), dwarf_all['apparent_magnitude_v'].max(), (dwarf_all['apparent_magnitude_v'].min(), dwarf_all['apparent_magnitude_v'].max()))
# dwarf_all = dwarf_all[(dwarf_all['apparent_magnitude_v'] >= apparent_magnitude_v[0]) & (dwarf_all['apparent_magnitude_v'] <= apparent_magnitude_v[1])]
# # filter by vlos systemic
# vlos_systemic = st.sidebar.slider('Vlos Systemic', dwarf_all['vlos_systemic'].min(), dwarf_all['vlos_systemic'].max(), (dwarf_all['vlos_systemic'].min(), dwarf_all['vlos_systemic'].max()))
# dwarf_all = dwarf_all[(dwarf_all['vlos_systemic'] >= vlos_systemic[0]) & (dwarf_all['vlos_systemic'] <= vlos_systemic[1])]


#-----create selections and conditions for interactivity-----#
color_scale = 'tableau10'
selection = alt.selection_point(fields=['source_pretty'], bind='legend',nearest=False,)

hover_selection = alt.selection_point(on='mouseover', nearest=False, empty=False)

color = alt.condition(
    hover_selection,
    alt.value('black'),
    alt.Color('source_pretty', scale=alt.Scale(scheme=color_scale), legend=alt.Legend(title='Source', symbolFillColor='black'), sort=table_names_pretty)
)
# sizeCondition=alt.condition(
#     hover_selection,
#     alt.SizeValue(100),
#     alt.SizeValue(0)
# )

strokeWidthCondition=alt.condition(
    hover_selection,
    alt.StrokeWidthValue(1),
    alt.StrokeWidthValue(0)
)

strokeErrorCondition=alt.condition(
    hover_selection,
    alt.value('black'),
    alt.Color('source_pretty', scale=alt.Scale(
            domain=table_names_pretty, scheme=color_scale), title='Source', legend=None)
)



#---------------------user select what values to show in tooltip----------------------#

with st.sidebar:
    def tooltip_items():
        tooltip_select = st.multiselect('What properties do you want to display in the tooltip?', valid_plot_cols, default=['name', 'host', plot_xaxis, plot_yaxis], format_func=lambda x: tab_desc[x]['label'])
        #tooltip = [alt.Tooltip(x, title=tab_desc[x]['label']) for x in tooltip_select]
        return tooltip_select
# print([tab_desc[x]['desc'] for x in tooltip_select])
#tooltip = tooltip_items()
    tooltip = [alt.Tooltip(x, title=tab_desc[x]['label']) for x in tooltip_items()]

    #st.write(tooltip)
# ---------------------plot---------------------- #


#print(tab_desc['ra'])
print(np.unique(master_df['source']))

# selection_x = alt.selection_interval(
#     bind='scales',
#     encodings=["x"],
#     zoom="wheel![event.altKey]",
# )
# selection_y = alt.selection_interval(
#     bind='scales',
#     encodings=["y"],
#     zoom="wheel![event.shiftKey]",
# )

# selection_both = alt.selection_interval(
#     bind='scales',
#     encodings=["x", "y"],
#     zoom="wheel![event.ctrlKey & !event.shiftKey & !event.altKey]",
# )

selection_x = alt.selection_interval(
    bind='scales',
    encodings=["x"],
    zoom="wheel![event.altKey]",
)

selection_y = alt.selection_interval(
    bind='scales',
    encodings=["y"],
    zoom="wheel![event.shiftKey]",
)

# selection_both = alt.selection_interval(
#     bind='scales',
#     encodings=["x", "y"],
#     zoom="wheel!",
# )

unique_sources = master_df['source_pretty'].unique()
unique_sources = sorted(unique_sources, key=lambda x: table_names_pretty.index(x))
print(unique_sources)

charts_to_layer = []
errors_to_layer = []

base_chart = alt.Chart(master_df).mark_point(filled=True, opacity=1, size=50).encode(
     x=alt.X(plot_xaxis, type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=xlabel), 
     y=alt.Y(plot_yaxis, type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=ylabel),
     color=alt.Color('source_pretty',scale=alt.Scale(scheme=color_scale, domain=table_names_pretty)),
                     
     tooltip = tooltip,
     #size=sizeCondition,
     strokeWidth=strokeWidthCondition,
     stroke=alt.value('black'),
     #size=alt.when(selection).then(alt.value(100)).otherwise(alt.value(0)),
     shape=alt.Shape('source_pretty', scale=alt.Scale(scheme=color_scale, domain=table_names_pretty)),
     #order=alt.Order('source_pretty'),
     ).add_params(hover_selection, selection_x, selection_y)#.transform_filter(selection)

#charts_to_layer.append(base_chart)

#help(alt.Chart.configure_point)

#print(plot_yaxis+"_ul")
if plot_yaxis+"_ul" in master_df.keys():
    tooltip.append(alt.Tooltip(plot_yaxis+"_ul", title=tab_desc[plot_yaxis+"_ul"]['label']))

    yul = alt.Chart(master_df).mark_point(shape="arrow", filled=True, strokeWidth=10).encode(
        x=alt.X(plot_xaxis, type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""),
        y=alt.Y(plot_yaxis+"_ul", type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
        #y2=alt.Y2(plot_yaxis+"_upper"),
        color=alt.Color('source_pretty', scale=alt.Scale(domain=table_names_pretty, scheme=color_scale), title='Source', legend=None),
        tooltip=tooltip,
        angle=alt.value(180),
        stroke=alt.value('black'),
        strokeWidth=strokeWidthCondition,
        #size=alt.value(50),
        yOffset=alt.value(10), # not sure why 10 works here to move the arrow to go from the center to the top of the point
        #size=alt.when(selection).then(alt.value(10)).otherwise(alt.value(0))
    ).transform_filter(
        (alt.datum[plot_yaxis+"_upper"] != 0) & (alt.datum[plot_xaxis] != 0) 
    )#.add_params(hover_selection)
    charts_to_layer.append(yul)

if plot_xaxis+"_ul" in master_df.keys():
    tooltip.append(alt.Tooltip(plot_xaxis+"_ul", title=tab_desc[plot_xaxis+"_ul"]['label']))

    xul = alt.Chart(master_df).mark_point(shape="arrow", filled=True).encode(
        x=alt.X(plot_xaxis+"_ul", type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""),
        y=alt.Y(plot_yaxis, type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
        size=alt.value(500),
        color=alt.Color('source_pretty', scale=alt.Scale(domain=table_names_pretty, scheme=color_scale), title='Source', legend=None),
        tooltip=tooltip,
        angle=alt.value(-90),
        xOffset=alt.value(-10),
        stroke=alt.value('black'),
        strokeWidth=strokeWidthCondition,
        #size=alt.when(selection).then(alt.value(10)).otherwise(alt.value(0))
    ).transform_filter(
        (alt.datum[plot_xaxis+"_upper"] != 0) & (alt.datum[plot_yaxis] != 0) 
    )#.add_params(hover_selection)
    charts_to_layer.append(xul)

if show_xerr:
    xerrorbars = alt.Chart(master_df).mark_errorbar(ticks=True).encode(
        x=alt.X(plot_xaxis+"_low", type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""),
        y=alt.Y(plot_yaxis, type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
        x2=alt.X2(plot_xaxis+"_high", title=""),
        size=alt.value(500),
        #stroke=alt.value('black'),
        #strokeWidth=strokeErrorCondition,
        color=strokeErrorCondition,
        # color=alt.Color('source_pretty', scale=alt.Scale(
        #     domain=table_names_pretty, scheme=color_scale), title='Source', legend=None),
        tooltip=alt.value(None),
        #opacity=alt.when(selection).then(alt.value(1)).otherwise(alt.value(0))
    ).transform_filter(
        (alt.datum[plot_xaxis+"_low"] != 0) & (alt.datum[plot_xaxis+"_high"] != 0) & (alt.datum[plot_yaxis] != 0)
    )#.add_params(hover_selection)
    charts_to_layer.append(xerrorbars)

    if plot_yaxis+"_ul" in master_df.keys():
        xup_errorbars = alt.Chart(master_df).mark_errorbar(ticks=True).encode(
            x=alt.X(plot_xaxis+"_low", type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""),
            y=alt.Y(plot_yaxis+"_ul", type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
            x2=alt.X2(plot_xaxis+"_high", title=""),
            size=alt.value(500),
            #stroke=alt.value('black'),
            color=strokeErrorCondition,
            #color=alt.Color('source_pretty', scale=alt.Scale(
            #domain=table_names_pretty, scheme=color_scale), title='Source', legend=None),
            tooltip=alt.value(None),
            #opacity=alt.when(selection).then(alt.value(1)).otherwise(alt.value(0))
        ).transform_filter(
            (alt.datum[plot_xaxis+"_low"] != 0) & (alt.datum[plot_xaxis+"_high"] != 0) & (alt.datum[plot_yaxis+"_ul"] != 'nan')
        )#.add_params(hover_selection)
        charts_to_layer.append(xup_errorbars)

if show_yerr:
    yerrorbars = alt.Chart(master_df).mark_errorbar(ticks=True).encode(
        x=alt.X(plot_xaxis, type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""), 
        y=alt.Y(plot_yaxis+"_low", type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
        y2=alt.Y2(plot_yaxis+"_high", title=""),
        # yError=alt.YError(plot_yaxis + '_low'),
        # yError2=alt.YError(plot_yaxis + '_high')
        #stroke=alt.value('black'),
        #strokeWidth=strokeWidthCondition,
        size=alt.value(500),
        color=strokeErrorCondition,
        # color=alt.Color('source_pretty', scale=alt.Scale(
        #     domain=table_names_pretty, scheme=color_scale), title='Source', legend=None),
        tooltip=alt.value(None),
        #opacity=alt.when(selection).then(alt.value(1)).otherwise(alt.value(0))
        ).transform_filter(
        (alt.datum[plot_yaxis+"_low"] != 0) & (alt.datum[plot_yaxis+"_high"] != 0) & (alt.datum[plot_xaxis] != 0)
        )#.add_params(hover_selection)
    charts_to_layer.append(yerrorbars)

    if plot_xaxis+"_ul" in master_df.keys():
        yup_errorbars = alt.Chart(master_df).mark_errorbar(ticks=True).encode(
            x=alt.X(plot_xaxis+"_ul", type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""), 
            y=alt.Y(plot_yaxis+"_low", type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
            y2=alt.Y2(plot_yaxis+"_high", title=""),
            size=alt.value(500),
            #stroke=alt.value('black'),
            #strokeWidth=strokeWidthCondition,
            color=strokeErrorCondition,
            # color=alt.Color('source_pretty', scale=alt.Scale(
            # domain=table_names_pretty, scheme=color_scale), title='Source', legend=None),
            tooltip=alt.value(None),
            #opacity=alt.when(selection).then(alt.value(1)).otherwise(alt.value(0))
            ).transform_filter(
            (alt.datum[plot_yaxis+"_low"] != 0) & (alt.datum[plot_yaxis+"_high"] != 0) & (alt.datum[plot_xaxis+"_ul"] != 'nan')
            )#.add_params(hover_selection)
        charts_to_layer.append(yup_errorbars)


charts_to_layer.append(base_chart)

# chart2 = alt.Chart(dsph_m31).mark_circle().encode(
#      x=alt.X(plot_xaxis, scale=alt.Scale(type=type_x, reverse=reverse_x), title=xlabel), 
#      y=alt.Y(plot_yaxis, scale=alt.Scale(type=type_y, reverse=reverse_y), title=tab_desc[plot_yaxis]['label']),
#      color=alt.Color('host', scale=alt.Scale(scheme='tableau20'), legend=alt.Legend(title='Host'))
#      ).interactive()
print(charts_to_layer)
#charts_to_layer = np.concatenate([charts_to_layer, errors_to_layer])
#print(charts_to_layer.reverse())
#st.altair_chart(errors_to_layer[0], use_container_width=True)
def plot_dwarf_all():
    #layered = base_chart+xerrorbars
    layered = alt.layer(*charts_to_layer).configure_legend(titleFontSize=18,labelFontSize=10).resolve_scale(shape='independent', color='independent').resolve_legend(color='independent', size='independent').configure_legend(symbolStrokeWidth=0)
    #st.altair_chart(alt.layer(*charts_to_layer).add_params(selection_x, selection_y).configure_legend(titleFontSize=18,labelFontSize=10).resolve_scale(shape='independent', color='independent').resolve_legend(color='independent', size='independent').configure_legend(symbolStrokeWidth=0), use_container_width=True)
    st.altair_chart(layered, use_container_width=True)
with st.container():
    plot_dwarf_all()

#st.altair_chart(charts_to_layer[1], use_container_width=True)

#st.toast('Welcome to the Local Volume Database!')





st.text('all dwarfs')
st.dataframe(master_df, use_container_width=True) # use_container_width doesn't work??
#st.text('dSphs in MW')
#st.dataframe(dsph_mw, use_container_width=False)
#st.dataframe(misc_host)

