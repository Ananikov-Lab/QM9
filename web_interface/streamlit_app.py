import streamlit as st

# Set white theme for Streamlit
st.set_page_config(page_title="Reaction Viewer")

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

# Custom CSS styles for modernized navigation menu
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

# Sidebar navigation with styled menu
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

# Page selection logic
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

