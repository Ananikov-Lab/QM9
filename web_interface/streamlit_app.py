import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
import base64

df_reactions = pd.read_csv('supporting_files/df_reactions.csv')

base_url = "static/dataset_png/"
df_reactions["Image"] = base_url + df_reactions["Image"].astype(str)

html_table = """
<html>
  <head>
    <!-- Подключаем CSS и JS DataTables с CDN -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.1/css/jquery.dataTables.css">
    <script src="https://code.jquery.com/jquery-3.5.1.js"></script>
    <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.13.1/js/jquery.dataTables.js"></script>
    <style>
      /* Общая стилистика */
      body {{
        background-color: white;
        color: black;
        font-family: 'Roboto', sans-serif;
        margin: 20px;
      }}
      
      /* Стилизация таблицы */
      table {{
        width: 100%;
        border-collapse: collapse;
        margin: 20px 0;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
      }}
      
      th, td {{
        padding: 12px 15px;
        border: 1px solid #ddd;
      }}
      
      th {{
        background-color: #f8f9fa;
        text-transform: uppercase;
        font-weight: bold;
      }}
      
      tr:hover {{
        background-color: #f1f1f1;
      }}
      
      img {{
        max-width: 150px;
        height: auto;
        border-radius: 8px;
      }}
      
      /* Стилизация элементов DataTables */
      .dataTables_wrapper .dataTables_filter input {{
          border: 1px solid #ddd;
          border-radius: 8px;
          padding: 8px;
      }}
      
      .dataTables_wrapper .dataTables_length select {{
          border: 1px solid #ddd;
          border-radius: 8px;
          padding: 8px;
      }}
    </style>
  </head>
  <body>
    <table id="myTable">
      <thead>
        <tr>
          <th>Reaction Image</th>
          <th>Energy (kcal/mol)</th>
        </tr>
      </thead>
      <tbody>
"""

def get_image_base64(image_path):
    with open(image_path, "rb") as img_file:
        b64_string = base64.b64encode(img_file.read()).decode("utf-8")
    return f"data:image/png;base64,{b64_string}"

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
        max-width: 900px;
        margin: auto;
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
    </ul>
    """,
    unsafe_allow_html=True
)


menu = st.query_params.get("nav", "Home")

if menu == "Home":
    st.title("Welcome to Reaction Viewer")
    st.markdown(
        """
        **The search for new chemical transformations is a key goal in modern chemistry.**
        
        The scope of available reactions defines the accessible chemical space and therefore can limit development in various molecule-dependent areas, such as drug discovery, medicine, agriculture, materials science, and life sciences, among others.
        
        This tool provides an interactive way to explore molecular transformations visually, helping in rapid assessment of chemical reactions and their thermodynamic feasibility.
        
        <h3>Authors:</h3>
        <p><strong>Nikita I. Kolomoets<sup>†,a</sup>, Daniil A. Boiko<sup>†,a</sup>, Leonid V. Romashov<sup>a</sup>, 
        Kirill S. Kozlov<sup>a</sup>, Evgeniy G. Gordeev<sup>a</sup>, Alexey S. Galushko<sup>a</sup>, 
        Valentine P. Ananikov<sup>a</sup>*</strong></p>

        <p><sup>a</sup> Zelinsky Institute of Organic Chemistry, Russian Academy of Sciences, Leninsky Prospekt 47, Moscow, 119991, Russia</p>
        <p><sup>†</sup> These authors contributed equally</p>
        <p><sup>*</sup> Corresponding author: val@ioc.ac.ru</p>
        
        <h3>How to use:</h3>
        <ul>
            <li>Navigate to the "Reaction Viewer" section.</li>
            <li>Select molecular display parameters.</li>
            <li>Analyze chemical reactions interactively.</li>
        </ul>
        """,
        unsafe_allow_html=True
    )
    
elif menu == "ReactionViewer":
    st.title("Reaction Viewer")
    st.markdown(
    """
    <iframe src="http://mc.zioc.su:8511" width="100%" height="1000px" frameborder="0"></iframe>
    """,
    unsafe_allow_html=True
    )
    
    st.write("### Reaction Table")

    rows_per_page = 100
    total_rows = len(df_reactions)
    total_pages = (total_rows - 1) // rows_per_page + 1


    page_number = st.number_input("Page Number", 
                                  min_value=1, 
                                  max_value=total_pages, 
                                  value=1, 
                                  step=1)


    start_idx = (page_number - 1) * rows_per_page
    end_idx = start_idx + rows_per_page

    df_page = df_reactions.iloc[start_idx:end_idx]


    for _, row in df_page.iterrows():
        image_data = get_image_base64(row['Image'])
        html_table += "<tr>"
        html_table += f"<td><img src='{image_data}' style='max-width:450px; height:auto;'/></td>"
        html_table += f"<td>{round(row['Energy (kcal/mol)'], 2)}</td>"
        html_table += f""
        html_table += "</tr>"

    html_table += """
        </tbody>
        </table>
        <script>
        $(document).ready(function () {
            $('#myTable').DataTable({
                "order": []  // начальное отсутствие сортировки
            });
        });
        </script>
    </body>
    </html>
    """

    components.html(html_table, height=600, scrolling=True)
