import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
import base64

df_reactions = pd.read_csv('supporting_files/df_reactions_multi.csv', header=[0,1])

base_url = "static/dataset_png/"
df_reactions[('Reaction', 'Image')] = base_url + df_reactions[('Reaction', 'Image')].astype(str)

def get_image_base64(image_path):
    with open(image_path, "rb") as img_file:
        b64_string = base64.b64encode(img_file.read()).decode("utf-8")
    return f"data:image/png;base64,{b64_string}"

html_table = """
<html>
  <head>
    <!-- Подключаем CSS и JS DataTables с CDN -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.1/css/jquery.dataTables.css">
    <script src="https://code.jquery.com/jquery-3.5.1.js"></script>
    <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.13.1/js/jquery.dataTables.js"></script>
    <style>
    [data-testid="stThemeSwitcher"] {
        display: none !important;
    }

    html, body {
        margin: 0;
        padding: 0;
        height: 100%;
        background-color: white;
        font-family: 'Roboto', sans-serif;
        color: black;
    }

    .table-container {
        width: 100% !important;
        margin: 0 !important;
        padding: 0 !important;
        box-sizing: border-box;
    }

    table {
        width: 100%;
        border-collapse: collapse;
        margin: 20px 0;
        table-layout: auto;
    }

    th, td {
        padding: 0.75rem;
        border-top: 1px solid #dee2e6;
    }

    th {
        background-color: #f8f9fa;
        color: #212529;
        text-transform: none;
        font-weight: 500;
        text-align: left;
        padding: 0.75rem;
        vertical-align: bottom;
        border-bottom: 2px solid #dee2e6;
    }

    td {
        border-left: 1px solid #dee2e6 !important;
        border-right: 1px solid #dee2e6 !important;
    }

    tr:hover {
        background-color: #f1f1f1;
        transition: background-color 0.3s ease;
    }

    tbody tr:nth-child(even) {
        background-color: #fafafa;
    }

    img {
        max-width: 150px;
        height: auto;
        border-radius: 8px;
    }

    .dataTables_wrapper .dataTables_filter input {
        border: 1px solid #ddd;
        border-radius: 8px;
        padding: 8px;
    }

    .dataTables_wrapper .dataTables_length select {
        border: 1px solid #ddd;
        border-radius: 8px;
        padding: 8px;
    }

    .block-container {
        padding-top: 0rem !important;
        margin-top: 0rem !important;
    }

    h1, h2, h3, h4, h5, h6 {
        margin-top: 0.2rem !important;
    }
    .table {
        width: 100%;
        margin-bottom: 1rem;
        color: #212529;
        border-collapse: separate;
    }
    .table thead th {
        vertical-align: bottom;
        border-bottom: 2px solid #dee2e6;
        font-size: 16px !important;
        padding: 4px !important;
    }
    .table-hover tbody tr:hover {
        background-color: #f1f1f1;
    }
    .dataTables_wrapper .dataTables_filter label {
        font-size: 0.85rem;
        color: #6c757d;
        display: inline-flex;
        align-items: center;
        gap: 5px;
    }
    .dataTables_wrapper .dataTables_filter label input {
        border-radius: 4px;
        border: 1px solid #ced4da;
        padding: 0.25rem 0.5rem;
    }

    .dataTables_scrollHeadInner, .dataTables_scrollBody table {
        width: 100% !important;
    }

    </style>
  </head>
  <body>
    <div class="table-container">
        <table id="myTable" class="table table-hover">

"""

html_table += "<thead>"
columns = df_reactions.columns
n_levels = columns.nlevels

for level in range(n_levels):
    html_table += "<tr>"
    col_idx = 0
    while col_idx < len(columns):
        current_col = columns[col_idx][level]
        if 'Unnamed' in current_col:
            current_col = ''

        colspan = 1
        for next_col in range(col_idx + 1, len(columns)):
            if columns[next_col][level] == columns[col_idx][level]:
                colspan += 1
            else:
                break

        rowspan = 1
        if colspan == 1:
            rowspan = 1
            for sub_level in range(level + 1, n_levels):
                if 'Unnamed' in columns[col_idx][sub_level]:
                    rowspan += 1
                else:
                    break

        rowspan_html = f" rowspan='{rowspan}'" if rowspan > 1 else ""
        colspan_html = f" colspan='{colspan}'" if colspan > 1 else ""

        html_table += f"<th scope='col'{rowspan_html}{colspan_html}>{current_col}</th>"
        col_idx += colspan
    html_table += "</tr>"
html_table += "</thead>"



st.set_page_config(page_title="Reaction Viewer", layout="wide")

st.markdown(
    """
    <style>
        iframe {
            display: block;
            border: none;
            width: 100%;
            height: 100vh;
        }
    </style>
    """,
    unsafe_allow_html=True
)

st.markdown(
    """
    <style>
    @import url("https://fonts.googleapis.com/css?family=Roboto:400,700");
    
    body {
        background-color: white;
        color: black;
        font-family: 'Roboto', sans-serif;
    }
    .reportview-container {
        background: white;
    }
    .stButton>button {
        background-color: #f0f0f0;
        color: black;
        border-radius: 8px;
        border: 1px solid #ddd;
        padding: 10px;
    }
    .stTextInput>div>div>input, .stSelectbox>div>div>div {
        border-radius: 8px;
        border: 1px solid #ddd;
        padding: 10px;
    }
    .block-container {
        margin: auto;
        padding-bottom: 0 !important;
    }
    .dataTables_wrapper .row {
        margin-bottom: 0 !important;
    }
    .sidebar .sidebar-content {
        background-color: #f8f9fa;
        border-right: 1px solid #ddd;
        padding: 20px;
        border-radius: 8px;
        width: 14.2%;
    }
    ul.nav-menu {
        list-style-type: none;
        padding: 0;
    }
    ul.nav-menu li {
        padding: 10px 0;
    }
    ul.nav-menu li a {
        display: block;
        padding: 8px 15px;
        font-size: 18px;
        font-weight: bold;
        text-decoration: none;
        text-transform: uppercase;
        color: black;
        border-radius: 6px;
        transition: background-color 0.3s ease-in-out;
    }
    ul.nav-menu li a:hover {
        background-color: black;
        color: white;
    }
    sup {
        vertical-align: super;
        font-size: smaller;
    }
    ul, p {
        text-align: justify;
        line-height: 1.3;
        margin-bottom: 5px;
    }
    h3 {
        text-align: center;
    }
    [data-testid="stAppViewContainer"],
    [data-testid="stSidebar"] {
        background-color: white !important;
    }
    .block-container {
        padding-top: 0rem !important;
        margin-top: 0rem !important;
    }

    h1, h2, h3, h4, h5, h6 {
        margin-top: 0.2rem !important;
    }
    .pagination-input input {
         border-radius: 4px !important;
         border: 1px solid #ced4da !important;
         padding: 0.375rem 0.75rem !important;
         text-align: center !important;
         font-size: 1rem !important;
         width: 80px !important;
         margin: auto;
         display: block;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.sidebar.title("Navigation")
st.sidebar.markdown(
    """
    <ul class="nav-menu">
        <li><a href="/" target="_self">Home</a></li>
        <li><a href="/?nav=ReactionViewer" target="_self">Reaction Viewer</a></li>
        <li><a href="/?nav=ReactionTable" target="_self">Reaction Table</a></li>
    </ul>
    """,
    unsafe_allow_html=True
)

menu = st.query_params.get("nav", "Home")

if menu == "Home":
    st.title("Welcome to Reaction Viewer")
    st.markdown(
        """
        <p>The search for new chemical transformations is a key goal in modern chemistry.</p>
        <p>The scope of available reactions defines the accessible chemical space and therefore can limit development in various molecule-dependent areas, such as drug discovery, medicine, agriculture, materials science, and life sciences, among others.</p>

        <p>This tool provides an interactive way to explore molecular transformations visually, helping in the rapid assessment of chemical reactions and their thermodynamic feasibility.</p>

        <h2>🔎 Looking for a specific reaction?</h2>
        <p>You can <strong>search for your reaction</strong> using the <strong>reaction space exploration tool</strong> and the <strong>reaction table</strong>.</p>
        <ul>
            <li>Use the <strong>search function</strong> to quickly locate reactions of interest.</li>
            <li>Browse through the <strong>reaction table</strong> to explore different transformations and their thermodynamic properties.</li>
            <li>Navigate to the <strong>"Reaction Viewer"</strong> section to analyze reactions interactively.</li>
        </ul>

        <h2>👩‍🔬 Authors</h2>
        <p><strong>Nikita I. Kolomoets<sup> †,a</sup>, Daniil A. Boiko<sup> †,a</sup>, Leonid V. Romashov<sup> a</sup>, Kirill S. Kozlov<sup> a</sup>, Evgeniy G. Gordeev<sup> a</sup>, Alexey S. Galushko<sup> a</sup>, Valentine P. Ananikov<sup> *,a</sup></strong></p>

        <p><sup>a</sup> Zelinsky Institute of Organic Chemistry, Russian Academy of Sciences, Leninsky Prospekt 47, Moscow, 119991, Russia</p>

        <p><sup>†</sup> These authors contributed equally</p>
        <p><sup>*</sup> Corresponding author: <a href="mailto:val@ioc.ac.ru">val@ioc.ac.ru</a></p>

        <h2>📌 How to use</h2>
        <ul>
            <li>Navigate to the <strong>"Reaction Viewer"</strong> section.</li>
            <li>Select <strong>molecular display parameters</strong>.</li>
            <li>Search and analyze <strong>chemical reactions interactively</strong>.</li>
        </ul>

        """,
        unsafe_allow_html=True
    )
    
elif menu == "ReactionViewer":
    st.title("Reaction Viewer")
    st.markdown(
        """
        <style>
          html, body { height: 100%; margin: 0; padding: 0; }
          .center-container { display: flex; justify-content: center; align-items: center; height: 100%; }
        </style>
        <div class="center-container">
          <iframe src="http://mc.zioc.su:8511" style="width:800px; height:800px; border:none;"></iframe>
        </div>
        """,
        unsafe_allow_html=True
    )

elif menu == "ReactionTable":
    st.title("Reaction Table")

    rows_per_page = 10
    total_rows = len(df_reactions)
    total_pages = (total_rows - 1) // rows_per_page + 1

    params = st.query_params
    page_param = params.get("page", None)
    if page_param is None:
        page_number = 1
    elif isinstance(page_param, list):
        try:
            page_number = int(page_param[0])
        except ValueError:
            page_number = 1
    else:
        try:
            page_number = int(page_param)
        except ValueError:
            page_number = 1

    prev_disabled = "disabled" if page_number <= 1 else ""
    next_disabled = "disabled" if page_number >= total_pages else ""
    prev_page = page_number - 1 if page_number > 1 else 1
    next_page = page_number + 1 if page_number < total_pages else total_pages

    pagination_links = f'<li class="page-item {"active" if page_number==1 else ""}"><a class="page-link" href="?nav=ReactionTable&page=1" target="_self">1</a></li>'

    if page_number > 3:
        pagination_links += '<li class="page-item disabled"><a class="page-link" href="#">...</a></li>'

    for p in [page_number - 1, page_number, page_number + 1]:
        if p > 1 and p < total_pages:
            if p == page_number:
                pagination_links += f'<li class="page-item active"><a class="page-link" href="?nav=ReactionTable&page={p}" target="_self">{p}</a></li>'
            else:
                pagination_links += f'<li class="page-item"><a class="page-link" href="?nav=ReactionTable&page={p}" target="_self">{p}</a></li>'

    if page_number < total_pages - 2:
        pagination_links += '<li class="page-item disabled"><a class="page-link" href="#">...</a></li>'

    if total_pages > 1:
        pagination_links += f'<li class="page-item {"active" if page_number==total_pages else ""}"><a class="page-link" href="?nav=ReactionTable&page={total_pages}" target="_self">{total_pages}</a></li>'

    pagination_html = f"""
    <style>
    .pagination-container {{
        margin-bottom: 20px;
    }}
    .pagination {{
        display: flex;
        padding-left: 0;
        list-style: none;
        border-radius: 0.25rem;
        margin: 0;
    }}
    .pagination li {{
        margin: 0 2px;
    }}
    .pagination li a {{
        background-color: #fff;
        color: #000;
        border: 1px solid #000;
        border-radius: 0.25rem;
        padding: 0.5rem 0.75rem;
        font-size: 0.85rem;
        text-decoration: none;
        transition: background-color 0.2s ease, color 0.2s ease;
    }}
    .pagination li a:hover {{
        background-color: #000;
        color: #fff;
    }}
    .pagination li.active a {{
        background-color: #000;
        color: #fff;
        border-color: #000;
    }}
    .pagination li.disabled a {{
        color: #6c757d;
        pointer-events: none;
        background-color: #fff;
        border-color: #dee2e6;
    }}
    </style>
    <div class="pagination-container">
      <div style="display: flex; justify-content: space-between; align-items: center;">
        <span style="font-size: 0.85rem; font-weight: 600;">Reaction Pages</span>
        <ul class="pagination justify-content-end">
          {pagination_links}
        </ul>
      </div>
    </div>
    """
    st.markdown(pagination_html, unsafe_allow_html=True)
    
    start_idx = (page_number - 1) * rows_per_page
    end_idx = start_idx + rows_per_page
    df_page = df_reactions.iloc[start_idx:end_idx]

    table_body = "<tbody>"
    for _, row in df_page.iterrows():
        table_body += "<tr>"
        for col in df_reactions.columns:
            if col == ('Reaction', 'Image'):
                img_src = get_image_base64(row[col])
                table_body += f"<td><img src='{img_src}' style='max-width:200px; height:auto;'/></td>"
            elif col == ('Reaction', 'Energy (kcal/mol)'):
                table_body += f"<td>{round(row[col], 2)}</td>"
            else:
                cell_value = row[col] if pd.notna(row[col]) else ""
                table_body += f"<td>{cell_value}</td>"
        table_body += "</tr>"
    table_body += "</tbody>"

    final_html = html_table + table_body + """
        </table>
        </div>
        <script>
        $(document).ready(function() {
        var table = $('#myTable').DataTable({
            "order": [],
            "paging": false,
            "scrollY": "calc(100vh - 250px)",
            "scrollCollapse": true,
            "searching": false
            });

            function fixTableWidth() {
                var head = $('.dataTables_scrollHeadInner table');
                var body = $('.dataTables_scrollBody table');
                head.css('width', body.width() + 'px');
            }

            fixTableWidth();
            $(window).resize(fixTableWidth);
        });
        </script>
    </body>
    </html>
    """

    components.html(final_html, height=900, scrolling=True)
