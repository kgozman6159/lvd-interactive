# ----------import necessary packages---------- #
import streamlit as st
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import altair as alt

from astropy.io import fits
import astropy.table as table
from astropy.io import ascii

from astropy import units as u 

import astropy.coordinates as coord

#import local_volume_database # idk what this does yet
#from collections import OrderedDict
st.set_page_config(layout="wide",
    page_title="Local Volume Database",
    page_icon="🌌",
    initial_sidebar_state="expanded",
    menu_items={
    'Report a bug': "https://github.com/kgozman6159/lvd-interactive/issues",
    'About': "Made with :heart: by Katya Gozman using Streamlit :streamlit:, python, Altair, and Pandas. Check out the data on [![Repo](https://badgen.net/badge/icon/GitHub?icon=github&label)](https://github.com/apace7/local_volume_database)!"
  }

)

import altair as alt
import re
from webcolors import name_to_hex

#from streamlit_extras.badges import badge 
from streamlit_theme import st_theme
#import streamlit.components.v1 as components

#rom vega_datasets import data



# cars = data.cars()
# chart = alt.Chart(cars).mark_circle().encode(
#         x=alt.X('Miles_per_Gallon', title=''),
#         y='Weight_in_lbs',
#         color='Origin'
# )

# text = alt.Chart().mark_text(
#     align="center",
#     baseline="top",
#     fontSize=11,
#     fontWeight=600,
#     #color='#007bff',
#     href='https://stackoverflow.com'
# ).encode(
#     x=alt.value(200),  # pixels from left
#     y=alt.value(322),  # pixels from top
#     text=alt.value("Miles_per_Gallon")
# )

# chart = text + chart
# chart['usermeta'] = {
#     "embedOptions": {
#         'loader': {'target': '_blank'}
#     }
# }


# st.altair_chart(chart, use_container_width=True)
#print("*@&#*@&3847@&*#&@#")
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
        
        if ("ref" in key):
            combined_df[key] = combined_df[key].fillna("No reference")
            # for val in combined_df[key]:
            #     print(val)
            #     match = re.search(r'\d+', val)
            #     if match:
            #         print(val[match.start():])
                

            combined_df["bibcode_"+key] = combined_df[key].apply(lambda x: "https://ui.adsabs.harvard.edu/abs/"+ x[re.search(r'\d+', x).start():] if re.search(r'\d+', x) else "N/A")

        combined_df['image'] = r"https://vega.github.io/vega-datasets/data/ffox.png"

        #combined_df = combined_df.iloc[::-1]
        # for key in dwarf_all.keys():
        #     if key.endswith('_em') or key.endswith('_ep'):
        #         dwarf_all[key].fillna(0, inplace=True)

    return dwarf_all, dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate, misc_host, combined_df

dwarf_all, dsph_mw, dsph_m31, dsph_lf, dsph_lf_distant, gc_ambiguous, gc_mw_new, gc_harris, gc_dwarf_hosted, gc_other, candidate, misc_host, master_df = load_data()

# ---------get info about master_df before filtering--------- #
TOTAL_NUM_SYSTEMS = len(master_df)
ALL_SYSTEM_NAMES = master_df['name'].unique()
ALL_SOURCES = master_df['source_pretty'].copy()
source_mapping = {name: idx for idx, name in enumerate(table_names_pretty)} # turn each system source (MW Dwarfs, M31 Dwarfs, etc. into a number 0-9)
all_source_indices = master_df['source_pretty'].copy().map(source_mapping) # for mapping to colors later

#ALL_SOURCES = master_df['source_pretty_num'].copy()
#print(len(dwarf_all),len(master_df), len(dsph_mw)+len(dsph_m31)+len(dsph_lf)+len(dsph_lf_distant)+len(gc_ambiguous)+len(gc_mw_new)+len(gc_harris)+len(gc_dwarf_hosted)+len(gc_other)+len(candidate),len(misc_host))
#st.dataframe(dwarf_all, use_container_width=True) # use_container_width doesn't work??
#st.dataframe(master_df, use_container_width=True)

theme = st_theme()
#print(theme)
#---------------------title and instructions----------------------#


st.markdown(
    """
<style>
div[data-testid="stDialog"] div[role="dialog"] {
    width: 80vw;
    
}
</style>
""",
    unsafe_allow_html=True,
)


# css="""
# <style>
#     [data-testid="stDialog"] div[role="dialog"]{
#         background: white;
#         /*background-image: url(https://cdn.esahubble.org/archives/images/newsfeature/heic1909a.jpg);
#         background-size: contain;
#         background-repeat: no-repeat;*/
#         # background-image: radial-gradient(circle, white,white,white,white,grey);
#         opacity: 1;
#         /*text-color: white;*/

#     }
# </style>
# """
# st.write(css, unsafe_allow_html=True)

@st.dialog("Welcome to the interactive Local Volume Database!", width='large')
def tutorial():
    string = """
    This website lets you plot different properties of dwarf galaxies and globular clusters in the Local Volume. 
    These properties are collated in the Local Volume Database, which has been compiled and is maintained by Andrew Pace."""

    st.header(string)
    
    col1, col2 = st.columns([1,1])
    with col1:
        st.markdown("""
        #### ⬅️ You can use the sidebar to adjust the plotting parameters.
        - :red-background[**x-axis, y-axis**] Select the property to plot on the x- and y-axes.
          - Hover over the :material/help: icon next to each property to see an explanation of what it is.
          - Some selected parameters will also let you display error bars for those quantities and/or change the scale of the axis from linear to logarithmic.
          - Arrow markers in the plot indicate that the value is an upper limit.
        - :red-background[**Source**] Filter the data by what type of systems you want to display.
        - :red-background[**Filters**] Filter the data by specific properties.
        - :red-background[**Tooltip**] Select the properties to display in the tooltip. If they are not in the data, they will be displayed as "null".
        - :red-background[**Color selection**] Change the color of each source in the plot.
                      
                    """)
    with col2:
        st.markdown("""
                
    #### :bar_chart: The main plot in the middle shows the properties of the objects in the database. 
                
    :three_button_mouse: Click and drag to :blue-background[**pan**] around the plot. 
                
    :mag_right: Use your mouse wheel and hold down the 'shift' and 'alt' ('option' on a Mac) keys simultaneously to :blue-background[**zoom**] in and out.
                You can also zoom only in the x or y directions by using the mouse wheel and holding down the 'alt' or 'shift' keys, respectively.
                    Note that zooming is only possible for quantitative axes. Double click anywhere on the plot to reset the zoom.

    :flying_saucer: :blue-background[**Hover**] over the points to see more information about each object. You can change the information displayed in the tooltip in the sidebar.
                    
     :pushpin: In the menu above the plot, use the dropdown to :blue-background[**select**] as many systems as you want by name to highlight in the plot. 
                    You can also type to search for a specific system. If the system does not have a value for the selected x and y axes,
                    it will not show up in the list.
                

    """)
        
    st.markdown("""
                Click on the :material/more_vert: icon in the top right of the sidebar to change appearance settings and report bugs.
                """)
        
    st.caption("""
    ### :bulb: If you ever want to read this information again, click on the "show tutorial" button in the sidebar. Thanks for visiting! :blush:
    """)

    

if "show_tutorial" not in st.session_state:
    st.session_state.show_tutorial = True

if st.session_state.show_tutorial:
    tutorial()
    #st.session_state.show_tutorial = False

if st.sidebar.button("Show Tutorial", on_click=lambda: st.session_state.update(show_tutorial=True)):
    st.session_state.show_tutorial = True
    

with st.sidebar:
    st.title("Parameters")
# ---------------------dictionary of column labels and descriptions---------------------- #
#print(dwarf_all['M_V_high'])
tab_desc = pd.read_csv('table_descriptions.csv', index_col='property', keep_default_na=False)
tab_desc['reference'] = tab_desc['reference'].replace('N/A', '')
#print(tab_desc['reference'])
#print(tab_desc.iloc[0].unit)
tab_desc = tab_desc.T.to_dict(index='property')
# print(tab_desc.keys())
#st.write(tab_desc)
# st.write(tab_desc['rhalf_physical']['desc'])
valid_plot_cols = ['ra', 'dec', 'name', 'host', 'confirmed_real', 
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
## M_V -> L_V
def lum(m_x, m_x_sun=4.83):
    return pow(10., -0.4*(m_x - m_x_sun) )

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

    if tab_desc[axis]['reference'] != "":
        ref = tab_desc[axis]['reference']
    else:
        ref = 'confirmed_real'

    return type_axis, reverse_axis, axis_label, channel, show_error, ref

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
        plot_xaxis = st.selectbox('x-axis', valid_plot_cols, key="my_key1", format_func=lambda x: tab_desc[x]['label'], label_visibility='visible', help=f"{tab_desc[st.session_state.my_key1]['desc']}", on_change=lambda: st.session_state.update(show_tutorial=False))
        #right.caption("---", help=tab_desc[plot_xaxis]['desc'])
        type_x, reverse_x, xlabel, channel_x, show_xerr, xref = get_axis_specs(plot_xaxis, 'xaxis')
    with st.container(border=True, key='ycont') as ycont:
        plot_yaxis = st.selectbox('y-axis', valid_plot_cols, key="my_key2", format_func=lambda x: tab_desc[x]['label'], help=f"{tab_desc[st.session_state.my_key2]['desc']}", on_change=lambda: st.session_state.update(show_tutorial=False))
        #right.caption("---", help=tab_desc[plot_yaxis]['desc'])
        type_y, reverse_y, ylabel, channel_y, show_yerr, yref = get_axis_specs(plot_yaxis, 'yaxis')
# st.selectbox(
#     "Make a selection",
#     ["Default", "Apple", "Banana", "Carrot"],
#     key="my_key",
#     help=f"You've selected {st.session_state.my_key}"
# )






# ---------------------filtering---------------------- #
#filter by source
source = st.sidebar.multiselect('Source', table_names_pretty, default=table_names_pretty, on_change=lambda: st.session_state.update(show_tutorial=False))
if source:
    master_df = master_df[master_df['source_pretty'].isin(source)]
print(source)

# ---------------------filtering by columns---------------------- #
with st.sidebar:
    with st.expander("Filtering"):
        st.markdown("### Filters")
        filter_container = st.container()
        filter_container.markdown("Add filters for different parameters:")
        filter_columns = st.multiselect('Select parameters to filter by', valid_plot_cols, format_func=lambda x: tab_desc[x]['label'], key="filter_columns", on_change=lambda: st.session_state.update(show_tutorial=False))

        
        filter_values = {}
        for col in filter_columns:
            if tab_desc[col]['dtype'] == 'float64' or tab_desc[col]['dtype'] == 'int64':
                min_val, max_val = master_df[col].min(), master_df[col].max()
                filter_values[col] = filter_container.slider(f'{tab_desc[col]["label"]} ({tab_desc[col]["unit"]})', min_val, max_val, (min_val, max_val), step=0.01, on_change=lambda: st.session_state.update(show_tutorial=False))
                #filter_values[col] = filter_container.select_slider(f'{tab_desc[col]["label"]} ({tab_desc[col]["unit"]})', sorted(master_df[col]), (min_val, max_val), on_change=lambda: st.session_state.update(show_tutorial=False))

            else:
                unique_vals = master_df[col].unique()
                filter_values[col] = filter_container.multiselect(f'{tab_desc[col]["label"]}', unique_vals, on_change=lambda: st.session_state.update(show_tutorial=False))
        #st.form_submit_button('Apply Filters', on_click=lambda: st.session_state.update(show_tutorial=False))


        for col, val in filter_values.items():
            if isinstance(val, tuple):
                master_df = master_df[(master_df[col] >= val[0]) & (master_df[col] <= val[1])]
            else:
                master_df = master_df[master_df[col].isin(val)]

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




#---------------------user select color for each source----------------------#

def get_luminance(hex_color):
    color = hex_color[1:]

    hex_red = int(color[0:2], base=16)
    hex_green = int(color[2:4], base=16)
    hex_blue = int(color[4:6], base=16)
    return hex_red * 0.2126 + hex_green * 0.7152 + hex_blue * 0.0722

hex_codes = ['#4c78a8', '#f58518', '#e45756', '#72b7b2', '#54a24b', '#FFC000', '#b279a2', '#ff9da6', '#9d755d', '#bab0ac']
with st.sidebar:
    with st.popover("Color selection", use_container_width=False, help='Change the color of each source in the plot'):
        cols = st.columns(2, vertical_alignment="top", gap='medium')
        range_ = [cols[x>4].color_picker('%s'%table_names_pretty[x],hex_codes[x], key="color%i"%x, on_change=lambda: st.session_state.update(show_tutorial=False)) for x in range(10)]




string = ""
for i in range(len(hex_codes)):
    luminance = get_luminance(hex_color=range_[i])

    if luminance < 140:
        col = "white"
    else:
        col = 'black'

    string+="""
    <style>
        span[data-baseweb="tag"][aria-label="%s, close by backspace"]{
            background-color: %s;
            color: %s;
        }
    </style>
    """%(table_names_pretty[i], range_[i], col)
st.markdown(string, unsafe_allow_html=True)

#---------------------user select what values to show in tooltip----------------------#

with st.sidebar:
    def tooltip_items():
        tooltip_select = st.multiselect('What properties do you want to display in the tooltip?', valid_plot_cols, default=['name', 'host', plot_xaxis, plot_yaxis], format_func=lambda x: tab_desc[x]['label'], on_change=lambda: st.session_state.update(show_tutorial=False))
        # if tooltip_select:
        #     st.toast(f"Selected {len(tooltip_select)} properties to display in the tooltip")
        #tooltip = [alt.Tooltip(x, title=tab_desc[x]['label']) for x in tooltip_select]
        return tooltip_select
# print([tab_desc[x]['desc'] for x in tooltip_select])
#tooltip = tooltip_items()
    tooltip = [alt.Tooltip(x, title=tab_desc[x]['label']) for x in tooltip_items()]

#---------------------user select and highlight a certain galaxy by name using st.multiselect----------------------#

def gal_search():

    non_nan_galaxies = master_df.dropna(subset=[plot_xaxis, plot_yaxis])['name'].tolist()
    #st.write("Galaxies with non-nan values for selected x and y axis columns:", non_nan_galaxies)
    highlight = st.multiselect('Highlight a galaxy by name', non_nan_galaxies, format_func=lambda x: x, on_change=lambda: st.session_state.update(show_tutorial=False))
    return highlight
selected_gals = gal_search()

filtered_df = master_df[master_df['name'].isin(selected_gals)]
    # print(highlight)
    # print(master_df['name'].unique())

#print(selected_gals)

string = ""
for i in range(TOTAL_NUM_SYSTEMS-1):
    source_index = all_source_indices[i]
    luminance = get_luminance(hex_color=range_[source_index])
    #print(luminance)

    if luminance < 140:
        col = "white"
    else:
        col = 'black'


    string+="""
    <style>
        span[data-baseweb="tag"][aria-label="%s, close by backspace"]{
            background-color: %s;
            color: %s
        }
    </style>
    """%(ALL_SYSTEM_NAMES[i], range_[source_index], col)
st.markdown(string, unsafe_allow_html=True)

# st.markdown("""
#     <style>
#         span[data-baseweb="tag"][aria-label="MW Dwarfs, close by backspace"]{
#             background-color: %s;
#         }
#         span[data-baseweb="tag"][aria-label="M31 Dwarfs, close by backspace"]{
#             background-color: %s;
#         }
#         span[data-baseweb="tag"][aria-label="option3, close by backspace"]{
#             background-color: %s;
#         }
#     </style>
#     """%(hex_codes[0], hex_codes[1], hex_codes[2]), unsafe_allow_html=True)


#-----create selections and conditions for interactivity-----#
#color_scale = st.selectbox('Color scale', ['viridis', 'inferno', 'plasma', 'magma', 'cividis', 'accent', 'category10', 'category20', 'category20b', 'category20c', 'dark2', 'paired', 'pastel1', 'pastel2', 'set1', 'set2', 'set3', 'tableau10', 'tableau20'])
selection = alt.selection_point(fields=['source_pretty'], bind='legend',nearest=False,)

#print(alt.Color(scale=alt.Scale(scheme='category10')).to_dict()['scale']['scheme'])
#default_colors = ['blue', 'orange', 'green', 'red', 'olive', 'brown', 'pink', 'darkgreen', 'purple', 'cyan']
#hex_codes = [name_to_hex(x) for x in default_colors]


#print(range_)
hover_selection = alt.selection_point(on='mouseover', nearest=False, empty=False)

# color = alt.condition(
#     hover_selection,
#     alt.value('black'),
#     alt.Color('source_pretty', scale=alt.Scale(scheme=color_scale), legend=alt.Legend(title='Source'), sort=table_names_pretty)
# )

# color = alt.when(hover_selection).then(alt.value('black')).otherwise(alt.Color('source_pretty:N', 
#                                                                                scale=alt.Scale(scheme=color_scale), 
#                                                                                legend=alt.Legend(title='Source'), 
#                                                                                sort=table_names_pretty))


color = alt.when(hover_selection).then(alt.value('black')).otherwise(alt.Color('source_pretty:N', 
                                                                               scale=alt.Scale(range=range_), 
                                                                               legend=alt.Legend(title='Source'), 
                                                                               sort=table_names_pretty))
# sizeCondition=alt.condition(
#     hover_selection,
#     alt.SizeValue(100),
#     alt.SizeValue(0)
# )

# strokeWidthCondition=alt.condition(
#     hover_selection,
#     alt.StrokeWidthValue(1),
#     alt.StrokeWidthValue(0)
# )
#print("THEME", theme)
if theme == None:
    strokeColor = alt.value('black')
    strokeErrorCondition=alt.when(hover_selection).then(strokeColor).otherwise(alt.Color('source_pretty:N', scale=alt.Scale(
            domain=table_names_pretty, range=range_), title='Source', legend=None))
elif theme['base'] == 'light':
    strokeColor = alt.value('black')
    strokeErrorCondition=alt.when(hover_selection).then(strokeColor).otherwise(alt.Color('source_pretty:N', scale=alt.Scale(
            domain=table_names_pretty, range=range_), title='Source', legend=None))
elif theme['base'] == 'dark':
    strokeColor = alt.value('white')
    strokeErrorCondition=alt.when(hover_selection).then(strokeColor).otherwise(alt.Color('source_pretty:N', scale=alt.Scale(
            domain=table_names_pretty, range=range_), title='Source', legend=None))
else:
    strokeColor = alt.value('black')
    strokeErrorCondition=alt.when(hover_selection).then(strokeColor).otherwise(alt.Color('source_pretty:N', scale=alt.Scale(
            domain=table_names_pretty, range=range_), title='Source', legend=None))
    
strokeWidthCondition=alt.when(hover_selection).then(alt.StrokeWidthValue(1)).otherwise(alt.StrokeWidthValue(0))

# strokeErrorCondition=alt.condition(
#     hover_selection,
#     alt.value('black'),
#     alt.Color('source_pretty', scale=alt.Scale(
#             domain=table_names_pretty, scheme=color_scale), title='Source', legend=None)
# )







# ---------------------plot---------------------- #


#print(tab_desc['ra'])
#print(np.unique(master_df['source']))

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

charts_to_layer = []
errors_to_layer = []

when_hover = alt.when(hover_selection, empty=False)
selection_click = alt.selection_point(empty=False, on='click', nearest=False)

if len(selected_gals)!=0:
    opacity = alt.when(
        alt.FieldOneOfPredicate(field='name', oneOf=selected_gals)).then(alt.value(1)).otherwise(alt.value(0.05))
else:
    opacity = alt.value(1)


base_chart = alt.Chart(master_df[::-1]).mark_point(filled=True, size=50).encode(
     x=alt.X(plot_xaxis, type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=xlabel), 
     y=alt.Y(plot_yaxis, type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=ylabel),
     color=alt.Color('source_pretty',scale=alt.Scale(domain=table_names_pretty, range=range_), legend=alt.Legend(title='System Type')),
     #opacity=alt.when(brush).then(alt.value(1)).otherwise(alt.value(0.05)),
     #opacity=opacity,
     tooltip = tooltip,
     #size=sizeCondition,
     strokeWidth=strokeWidthCondition,
     stroke=strokeColor,
     #order=alt.value(0),
     #href=alt.when(alt.FieldOneOfPredicate(field=xref, oneOf=[0,1])).then(alt.Href("www.google.com:N")).otherwise(xref),
     #size=alt.when(selection).then(alt.value(100)).otherwise(alt.value(0)),
     shape=alt.Shape('source_pretty', scale=alt.Scale(domain=table_names_pretty), legend=alt.Legend(title='System Type')),
     #order=alt.Order('source_pretty'),
     ).add_params(selection_click, hover_selection, selection_x, selection_y)#.transform_filter(selection)

if plot_xaxis == 'M_V' and plot_yaxis == 'metallicity':
    x = np.arange( -20,4, .1)
    mass_met_data = pd.DataFrame({
        "x":x,
        'f(x)':-1.68 + 0.29 * np.log10(lum(x)/1e6)
    })

    mass_met = alt.Chart(mass_met_data).mark_line().encode(
        x='x',
        y='f(x)',
        color=alt.value("##FFAA00"),
        text=alt.value('Mass-Met'),
        tooltip=alt.value('Relation from Simon 2019'),

        
    )


    charts_to_layer.append(mass_met)


if plot_xaxis == 'rhalf_sph_physical' and plot_yaxis == 'M_V':
    x = np.arange(1e0, 1e4, 1)
    x = np.logspace(0, 4, 1000)
    for mu in [24, 26, 28, 30, 32]:
        const_mu = pd.DataFrame({
        "x":(x),
        'f(x)': mu - 36.57 - 2.5 * np.log10(2.*np.pi*(x/1000)**2)
        })
        #plt.plot(x,  [const_mu(mu, i/1000.) for i in x], c='k', lw=2, ls=':')
        const_mu_chart = alt.Chart(const_mu).mark_line().encode(
            x=alt.X('x',scale=alt.Scale(type=type_x, reverse=reverse_x)),
            y=alt.Y('f(x)', scale=alt.Scale(type=type_y, reverse=reverse_y)),
            color=alt.value("black"),
            text=alt.value(f'μᵥ = {mu} mag/arcsec²'),
            tooltip=alt.value(f'μᵥ= {mu} mag/arcsec²'),
            
        )
        charts_to_layer.append(const_mu_chart)



#charts_to_layer.append(base_chart)

#help(alt.Chart.configure_point)

#print(plot_yaxis+"_ul")

# text = alt.Chart(master_df).mark_text(
#     align="left", baseline="top", href='www.google.com',
# ).encode(
#     x=alt.value(0.0),
#     y=alt.value(3),
#     #x=alt.value(5),  # pixels from left
#     #y=alt.value(5),  # pixels from top
#     text=xref,
   
#     opacity=alt.condition(selection_click, alt.value(1.0), alt.value(0.0)),
#     size=alt.value(10),
#     )
#charts_to_layer.append(text)

#print(master_df['image'][0])
# image = alt.Chart(master_df).mark_image(
#     width=50, height=50).encode(
#         x=alt.value(0),
#         y=alt.value(3),
#         url="img:N", 
#         size=alt.value(1),
#         #opacity=alt.condition(selection_click, alt.value(1.0), alt.value(0.0)),
#         #href='xref'
#         )
#charts_to_layer.append(image)
# source = pd.DataFrame.from_records(
#     [
#         {
#             "x": 0.5,
#             "y": 0.5,
#             "img": "https://vega.github.io/vega-datasets/data/ffox.png",
#         },
#         {
#             "x": 1.5,
#             "y": 1.5,
#             "img": "https://vega.github.io/vega-datasets/data/gimp.png",
#         },
#         {
#             "x": 2.5,
#             "y": 2.5,
#             "img": "https://vega.github.io/vega-datasets/data/7zip.png",
#         },
#     ]
# )

# image = alt.Chart(source).mark_image(width=50, height=50).encode(x=alt.value(0.0), y=alt.value(3), url="img")

#st.altair_chart(text, use_container_width=True)

if show_xerr:
    xerrorbars = alt.Chart(master_df[::-1]).mark_errorbar(ticks=True).encode(
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
        order=alt.value(1)
        #opacity=alt.when(selection).then(alt.value(1)).otherwise(alt.value(0))
    ).transform_filter(
        (alt.datum[plot_xaxis+"_low"] != 0) & (alt.datum[plot_xaxis+"_high"] != 0) & (alt.datum[plot_yaxis] != 0)
    )#.add_params(hover_selection)
    charts_to_layer.append(xerrorbars)

    if plot_yaxis+"_ul" in master_df.keys():
        xup_errorbars = alt.Chart(master_df[::-1]).mark_errorbar(ticks=True).encode(
            x=alt.X(plot_xaxis+"_low", type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""),
            y=alt.Y(plot_yaxis+"_ul", type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
            x2=alt.X2(plot_xaxis+"_high", title=""),
            size=alt.value(500),
            #stroke=alt.value('black'),
            color=strokeErrorCondition,
            #color=alt.Color('source_pretty', scale=alt.Scale(
            #domain=table_names_pretty, scheme=color_scale), title='Source', legend=None),
            tooltip=alt.value(None),
            order=alt.value(1)
            #opacity=alt.when(selection).then(alt.value(1)).otherwise(alt.value(0))
        ).transform_filter(
            (alt.datum[plot_xaxis+"_low"] != 0) & (alt.datum[plot_xaxis+"_high"] != 0) & (alt.datum[plot_yaxis+"_ul"] != 'nan')
        )#.add_params(hover_selection)
        charts_to_layer.append(xup_errorbars)

if show_yerr:
    yerrorbars = alt.Chart(master_df[::-1]).mark_errorbar(ticks=True).encode(
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
        order=alt.value(1)
        #opacity=alt.when(selection).then(alt.value(1)).otherwise(alt.value(0))
        ).transform_filter(
        (alt.datum[plot_yaxis+"_low"] != 0) & (alt.datum[plot_yaxis+"_high"] != 0) & (alt.datum[plot_xaxis] != 0)
        )#.add_params(hover_selection)
    charts_to_layer.append(yerrorbars)

    if plot_xaxis+"_ul" in master_df.keys():
        yup_errorbars = alt.Chart(master_df[::-1]).mark_errorbar(ticks=True).encode(
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
            order=alt.value(1)
            #opacity=alt.when(selection).then(alt.value(1)).otherwise(alt.value(0))
            ).transform_filter(
            (alt.datum[plot_yaxis+"_low"] != 0) & (alt.datum[plot_yaxis+"_high"] != 0) & (alt.datum[plot_xaxis+"_ul"] != 'nan')
            )#.add_params(hover_selection)
        charts_to_layer.append(yup_errorbars)

if plot_yaxis+"_ul" in master_df.keys():
    tooltip.append(alt.Tooltip(plot_yaxis+"_ul", title=tab_desc[plot_yaxis+"_ul"]['label']))

    yul = alt.Chart(master_df[::-1]).mark_point(shape="arrow", filled=True, strokeWidth=10).encode(
        x=alt.X(plot_xaxis, type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""),
        y=alt.Y(plot_yaxis+"_ul", type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
        #y2=alt.Y2(plot_yaxis+"_upper"),
        color=alt.Color('source_pretty', scale=alt.Scale(domain=table_names_pretty, range=range_), title='Source', legend=None),
        tooltip=tooltip,
        angle=alt.value(180),
        stroke=strokeColor,
        strokeWidth=strokeWidthCondition,
        size=alt.value(500),
        order=when_hover.then(alt.value(1)).otherwise(alt.value(0)),
        yOffset=alt.value(10), # not sure why 10 works here to move the arrow to go from the center to the top of the point
        #size=alt.when(selection).then(alt.value(10)).otherwise(alt.value(0))
    ).transform_filter(
        (alt.datum[plot_yaxis+"_upper"] != 0) & (alt.datum[plot_xaxis] != 0) 
    )#.add_params(hover_selection)
    charts_to_layer.append(yul)

if plot_xaxis+"_ul" in master_df.keys():
    tooltip.append(alt.Tooltip(plot_xaxis+"_ul", title=tab_desc[plot_xaxis+"_ul"]['label']))

    xul = alt.Chart(master_df[::-1]).mark_point(shape="arrow", filled=True).encode(
        x=alt.X(plot_xaxis+"_ul", type=channel_x, scale=alt.Scale(type=type_x, reverse=reverse_x), title=""),
        y=alt.Y(plot_yaxis, type=channel_y, scale=alt.Scale(type=type_y, reverse=reverse_y), title=""),
        size=alt.value(500),
        color=alt.Color('source_pretty', scale=alt.Scale(domain=table_names_pretty, range=range_), title='Source', legend=None),
        tooltip=tooltip,
        angle=alt.value(-90),
        xOffset=alt.value(-10),
        stroke=strokeColor,
        strokeWidth=strokeWidthCondition,
        #size=alt.when(selection).then(alt.value(10)).otherwise(alt.value(0))
    ).transform_filter(
        (alt.datum[plot_xaxis+"_upper"] != 0) & (alt.datum[plot_yaxis] != 0) 
    )#.add_params(hover_selection)
    charts_to_layer.append(xul)


charts_to_layer.append(base_chart)
#charts_to_layer.append(filter_chart)


# chart2 = alt.Chart(dsph_m31).mark_circle().encode(
#      x=alt.X(plot_xaxis, scale=alt.Scale(type=type_x, reverse=reverse_x), title=xlabel), 
#      y=alt.Y(plot_yaxis, scale=alt.Scale(type=type_y, reverse=reverse_y), title=tab_desc[plot_yaxis]['label']),
#      color=alt.Color('host', scale=alt.Scale(scheme='tableau20'), legend=alt.Legend(title='Host'))
#      ).interactive()
#print(charts_to_layer)
#charts_to_layer = np.concatenate([charts_to_layer, errors_to_layer])
#print(charts_to_layer.reverse())
#st.altair_chart(errors_to_layer[0], use_container_width=True)

def plot_dwarf_all():
    #layered = base_chart+xerrorbars
    layered = (alt.layer(*charts_to_layer).encode(opacity=opacity, order=when_hover.then(alt.value(1)).otherwise(alt.value(0)))).configure_legend(titleFontSize=18,labelFontSize=10).resolve_scale(shape='independent', color='independent').resolve_legend(color='independent', size='independent').configure_legend(symbolStrokeWidth=0)
    #print(layered.to_json())
    #print("@#&@*(#@(#)@)*#@#@")
    
    #st.altair_chart(alt.layer(*charts_to_layer).add_params(selection_x, selection_y).configure_legend(titleFontSize=18,labelFontSize=10).resolve_scale(shape='independent', color='independent').resolve_legend(color='independent', size='independent').configure_legend(symbolStrokeWidth=0), use_container_width=True)
    st.altair_chart(layered, use_container_width=True)#, on_select="rerun")
    #print(layered.to_dict()['params'])
    #st.write(event)
    #return event
with st.container():
    plot_dwarf_all()

st.markdown("Filtered Galaxies")
st.dataframe(filtered_df, use_container_width=True, selection_mode='multi-row', hide_index=False, on_select="rerun")

with st.expander("Description of Catalogs"):
    #st.markdown("### Description of Catalogs")
    st.markdown("""
            These descriptions are taken from the official [documentation on using Local Volume Database catalogs](https://local-volume-database.readthedocs.io/en/latest/usage.html#decription-of-catalogs-tables).
            - **MW Dwarfs:** Milky Way dwarf galaxies (the most distant dwarf galaxy is Eridanus II at ~ 350 kpc).
            - **M31 Dwarfs:** M31 dwarf galaxies
            - **Local Field Dwarfs:** dwarf galaxies outside of MW/M31 within the Local Field to a distance of ~ 3 Mpc. This is an extension of galaxies from McConnachie 2012 compilation.
            - **Local Field Distance Dwarfs:** dwarf galaxies in the Local Volume with distance > 3 Mpc. The limiting distance is set to ~10-40 Mpc (the approximate limits of HST/JWST). This table is not complete to known systems (it is complete for known systems to a distance < 3.5 Mpc).
            - **Ambiguous GCs:** systems with an ambiguous classification (referred to as ambiguous or hyper-faint compact stellar systems in the LVDB). These are all MW halo systems.
            - **New MW GCs:** newly discovered globular clusters or candidate globular clusters (i.e. post-Harris catalog). Many systems are at low Galactic latitude (abs(b) <10-20 deg) and candidate systems may be open clusters.
            - **Harris GCs:** globular clusters in Harris catalog (this excludes Koposov 1 and 2 which are in the gc_abmiguous table).
            - **Dwarf Hosted GCs:** globular clusters hosted by dwarf galaxies. This does not include the Sagittarius globular clusters which are in gc_harris. This catalog is incomplete for known systems.
            - **Other GCs:** for other globular clusters. (mostly for future work)
            - **Candidates:** known false-positive candidates, background galaxies, or low confidence candidates. Only included in the release page.
                        """)



    
st.markdown("Check out the data on [![Repo](https://badgen.net/badge/icon/GitHub?icon=github&label)](https://github.com/apace7/local_volume_database)!")
st.markdown("Look at the documentation on [![readthedocs](https://img.shields.io/badge/readthedocs-ffffff?logo=readthedocs&style=flat&color=ffffff&logoColor=8CA1AF)](https://local-volume-database.readthedocs.io/en/latest/index.html)")
st.markdown("Read the paper on [![arXiv](https://img.shields.io/badge/arXiv-ffffff?logo=arxiv&style=flat&color=ffffff&logoColor=B31B1B)](https://arxiv.org/abs/2411.07424)")
#st.altair_chart(charts_to_layer[1], use_container_width=True)


#badge(type="github", name="apace7/local_volume_database/")
# Add Link to your repo
# '''
#     Check out the data! [![Repo](https://badgen.net/badge/icon/GitHub?icon=github&label)](https://github.com/apace7/local_volume_database/) 

# '''
#st.markdown("Check out the data! [![Repo](https://badgen.net/badge/icon/GitHub?icon=github&label)](https://github.com/apace7/local_volume_database/) ",unsafe_allow_html=True)


# st.text('all dwarfs')
# st.dataframe(master_df, use_container_width=True) # use_container_width doesn't work??
#st.text('dSphs in MW')
#st.dataframe(dsph_mw, use_container_width=False)
#st.dataframe(misc_host)

