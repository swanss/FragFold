import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

'''
The following functions plot generic fragment-specific values as a function of position in the sequence
'''

def plotRawValuesOnSingle(df:pd.DataFrame,
                  x='fragment_center_aa',
                  y='weighted_contacts',
                  hue='protein_name',
                  max_only=False):
    '''Plot fragment values together on the same axes
    '''
    if max_only:
        df = df.iloc[df.groupby([hue,x])[y].idxmax()]
    aspect_ratio = (2 / 350) * max(150,(df[x].max() - df[x].min())) #based on FtsZ
    plt.figure(figsize=(5*aspect_ratio,5))
    ax = sns.lineplot(data=df,x=x,y=y,hue=hue)
    ymax = max(10,df[y].max())+5
    ax.set_ylim(0,ymax)
    ax.set_title('Average weighted contacts vs. fragment position')
    ax.set_ylabel('Average weighted contacts')
    ax.set_xlabel('Fragment position in protein (center amino acid)')
    return ax
    
def plotRawValuesOnFacetGrid(df:pd.DataFrame,
                  x='fragment_center_aa',
                  y='weighted_contacts',
                  row='protein_name',
                  max_only=False):
    '''Plot fragment values separately, but share the x-axis
    '''
    if df['fragment_parent_name'].nunique() != 1:
        raise ValueError("The fragments should be derived from a single protein")
    fragment_parent_name = df['fragment_parent_name'].iloc[0]
    if max_only:
        df = df.iloc[df.groupby([row,x])[y].idxmax()]
    aspect_ratio = (3 / 350) * max(150,(df[x].max() - df[x].min())) #based on FtsZ
    g = sns.relplot(data=df,
            x=x,y=y,row=row,
            kind='line',
            height=6,aspect=aspect_ratio,
            facet_kws={'sharey': False, 'sharex': True})
    conditions = df[row].unique()
    for i in range(g.axes.shape[0]):
        condition = conditions[i]
        assert f'{row} = ' + condition == g.axes[i,0].title.get_text()
        ymax = max(10,df[df[row]==condition][y].max())+5
        g.axes[i,0].set_ylim(0,ymax)
    g.figure.suptitle(f'Fragments from {fragment_parent_name} vs different protein partners')
    g.figure.subplots_adjust(top=.9)
    g.set_titles('Protein = {row_name}')
    g.set_ylabels('Average weighted contacts')
    g.set_xlabels('Fragment position in protein (center amino acid)')
    return g
