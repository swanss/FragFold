import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

'''
The following functions plot generic fragment-specific values as a function of position in the sequence
'''

def plotRawValuesOnSingle(df:pd.DataFrame,
                  x='fragment center (aa)',
                  y='weighted_contacts',
                  hue='condition',
                  max_only=False):
    '''Plot fragment values together on the same axes
    '''
    if max_only:
        df = df.iloc[df.groupby([hue,x])[y].idxmax()]
    aspect_ratio = (2 / 350) * (df[x].max() - df[x].min()) #based on FtsZ
    plt.figure(figsize=(5*aspect_ratio,5))
    ax = sns.lineplot(data=df,x='fragment center (aa)',y='weighted_contacts',hue='condition')
    ymax = max(10,df[y].max())+5
    ax.set_ylim(0,ymax)
    return ax
    
def plotRawValuesOnFacetGrid(df:pd.DataFrame,
                  x='fragment center (aa)',
                  y='weighted_contacts',
                  row='condition',
                  max_only=False):
    '''Plot fragment values separately, but share the x-axis
    '''
    if max_only:
        df = df.iloc[df.groupby([row,x])[y].idxmax()]
    aspect_ratio = (3 / 350) * (df[x].max() - df[x].min()) #based on FtsZ
    g = sns.relplot(data=df,
            x=x,y=y,row=row,
            kind='line',
            height=2,aspect=aspect_ratio,
            facet_kws={'sharey': False, 'sharex': True})
    conditions = df[row].unique()
    for i in range(g.axes.shape[0]):
        condition = conditions[i]
        assert 'condition = ' + condition == g.axes[i,0].title.get_text()
        ymax = max(10,df[df[row]==condition][y].max())+5
        g.axes[i,0].set_ylim(0,ymax)
    return g
