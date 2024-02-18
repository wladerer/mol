
import pandas as pd
import plotly.express as px

def plot_uv_vis_from_df(df: pd.DataFrame, color_key: str | None = None):
    fig = px.line(df, x="Excitation Energy (eV)", y="Intensity", markers=True, color=color_key)
    fig.show()
    


