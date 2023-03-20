import pandas as pd
from typing import List, Callable
import matplotlib.pyplot as plt

def ideal_gas_law(gas_constant: float, list_of_volumes: List, list_of_temperatures: List) -> pd.DataFrame:
    """_summary_

    Args:
        gas_constant (float): _description_
        list_of_volumes (List): _description_
        list_of_temperatures (List): _description_

    Returns:
        _type_: _description_
    """
    df = pd.DataFrame(columns=list_of_temperatures, index=list_of_volumes)
    df = df.apply(lambda x: (float(x.name) * gas_constant) / x.index)
    return df


# Plot Function for Task 1 - Ideal Gas Law

def plot_ideal_gas_law(ideal_gas_law_df: pd.DataFrame) -> None:
    for column in ideal_gas_law_df:
        plt.plot(ideal_gas_law_df.index, ideal_gas_law_df[column], label=column)
        

# Plot Function for Task 3 - Ideal Gas Law

def peng_robinson_eos(gas_constant: float, list_of_volumes: List, list_of_temperatures: List, bc: float, ac: float, alpha: Callable) -> pd.DataFrame:
    df = pd.DataFrame(columns=list_of_temperatures, index=list_of_volumes)
    temperature_aT = {T: (ac * alpha(T)) for T in list_of_temperatures}
    df = df.apply(lambda x: ((gas_constant * float(x.name)) / (x.index - bc)) - (temperature_aT[float(x.name)] / ((x.index * (x.index + bc) + bc * (x.index - bc)))), axis=0)
    return df