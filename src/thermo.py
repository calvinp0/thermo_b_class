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


# Plot Function for Task 2 - Ideal Gas Law

def plot_ideal_gas_law(ideal_gas_law_df: pd.DataFrame) -> None:
    for column in ideal_gas_law_df:
        plt.plot(ideal_gas_law_df.index, ideal_gas_law_df[column], label=column)
    plt.ylim(0, 15)
        

# Peng-Robinson EOS Function for Task 3

def peng_robinson_eos(gas_constant: float, list_of_volumes: List, list_of_temperatures: List, bc: float, ac: float, alpha: Callable) -> pd.DataFrame:
    df = pd.DataFrame(columns=list_of_temperatures, index=list_of_volumes)
    temperature_aT = {T: (ac * alpha(T)) for T in list_of_temperatures}
    df = df.apply(lambda x: ((gas_constant * float(x.name)) / (x.index - bc)) - (temperature_aT[float(x.name)] / ((x.index * (x.index + bc) + bc * (x.index - bc)))), axis=0)
    return df

# Plot Function for Task 4 - Peng-Robinson EOS

def plot_peng_robinson_eos(peng_robinson_eos_df: pd.DataFrame) -> None:
    for column in peng_robinson_eos_df:
        plt.plot(peng_robinson_eos_df.index, peng_robinson_eos_df[column], label=column)

# PLot Function for Task 6

def plot_example_ideal_gas_law_and_peng_robinson_eos() -> None:
    example = pd.read_csv('src/preos_df_example.csv', index_col=0)
    plot_peng_robinson_eos(example)
    plt.axhline(y=4.6, color='red', linestyle='--', label='$P^{0}$ at 284.0 K (Lit)')
    plt.xlabel('Volume (m^3)')
    plt.ylabel('Pressure (MPa)')
    plt.ylim(0, 15)
    plt.title('Peng-Robinson EOS - CO2 Example')
    plt.legend()
    plt.show()


def get_chemical_values(chemical_df: pd.DataFrame, molecules: List, column: str) -> List:
    if column not in chemical_df.columns:
        raise ValueError('Column not in DataFrame. Did you input the correct column name?')
    elif molecules is None:
        return ValueError('No molecules were inputted. Please input a list of molecules.')
    elif not isinstance(molecules, list):
        raise TypeError('Molecules must be a list.')
    elif [molecule in chemical_df['Element'].unique() for molecule in molecules] != [True, True]:
        raise ValueError('Molecules not in DataFrame. Did you input the correct molecules?')
    elif chemical_df is None:
        raise ValueError('No DataFrame was inputted. Please input the required DataFrame.')
    else:
        return chemical_df.loc[chemical_df['Element'].isin(molecules), column].values.tolist()
    

# Calculate ai for Task 7.1.1
def calculate_ai(chemical_df: pd.DataFrame, molecules: List, temperature: float) -> List:
    """_summary_

    Args:
        chemical_df (pd.DataFrame): _description_
        molecules (List): _description_
        temperature (float): _description_

    Returns:
        List: _description_
    """
    if temperature <= 0:
        raise ValueError('Temperature must be greater than 0.')
    elif not isinstance(temperature, float):
        raise TypeError('Temperature must be a float.')
    else:
        #  0.45724 * R**2 * df.loc['Tc_K']**2 / df.loc['Pc_MPa'] * (1 + df.loc['ki'] * (1 - (T/df.loc['Tc_K'])**0.5))**2
        return 0.45724 * (8.3145**2) * chemical_df.loc['Tc_K']**2/chemical_df.loc['Pc_MPa'] * (1 + chemical_df.loc['ki'] * (1 - (temperature/chemical_df.loc['Tc_K'] )**0.5))**2


# Task 7.2

def calculate_df_aij(molecules: List, df: pd.DataFrame, df_interaction: pd.DataFrame, df_aij: pd.DataFrame):
    for i in molecules:
        for j in molecules:
            df_aij.loc[i, j] = (df.loc['ai', i] * df.loc['ai', j])**0.5 * (1 - df_interaction.loc[i, j])
    return df_aij