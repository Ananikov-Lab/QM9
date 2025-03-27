import molplotly as mp
import plotly.express as px
import pandas as pd
import pickle

with open('supporting_files/dict_variables_qm9_1310.pkl', "rb") as f:
    data_dict = pickle.load(f)

df_final = data_dict


x = df_final["visual_tSNE_1024"][:, 0]
y = df_final["visual_tSNE_1024"][:, 1]

smiles_cols = ["smiles_reag", "smiles_prod"]

df_vis = pd.DataFrame({
    "x": x,
    "y": y,
    "smiles_reag": df_final["smiles_reag"],
    "smiles_prod": df_final["smiles_prod"],
})

fig_scatter = px.scatter(df_vis,
                         x="x",
                         y="y",
                         title="",
                         width=800,
                         height=800,
                         )
                         

fig_scatter = mp.add_molecules(fig=fig_scatter, 
                               df=df_vis, 
                               smiles_col=smiles_cols,
                               show_coords=False,
                               fontfamily='Roboto')
if __name__ == "__main__":
    fig_scatter.run_server(debug=True, host="0.0.0.0", port=5011)
