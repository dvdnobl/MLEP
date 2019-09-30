import pandas as pd
import scipy as sp

#reading datasets as csv files into panda dataframes
RNA_binding = pd.read_csv('priorities.csv')
fragment_peaks = pd.read_csv('joint-frag-mle-peak.csv')

#merging datasets bv systematic gene name, separating by screen 
merged = pd.merge(RNA_binding, fragment_peaks, on='yorf', how='inner')
merged.drop(['beckmannHi','beckmannLo','hogan','tsvetanova','mitchell','RNAbinding','essential','desc_x','desc_y','mlePeak','gene_y'], axis=1, inplace=True)
screen_15 = merged.drop(['mlePeak.18','nread.18','nbc.18'], axis=1)
screen_18 = merged.drop(['mlePeak.15','nread.15','nbc.15'], axis=1)

#taking the absolute value of the MLEP scores 
screen_15['mlePeak.15'] = screen_15['mlePeak.15'].apply(abs)   
screen_18['mlePeak.18'] = screen_18['mlePeak.18'].apply(abs)  

#MLEP from screen 15 into RNA-binding (score >= 2) and non-binding (score <= 1)
dist_15_binding = screen_15.loc[screen_15['score'] >= 2]
dist_15_nonbinding = screen_15.loc[screen_15['score'] <= 1]

#KDE plots for screen 15 
binding_plot_15 = dist_15_binding['mlePeak.15'].plot.kde()
nonbinding_plot_15 = dist_15_nonbinding['mlePeak.15'].plot.kde()

#To see distributions for screen 18, uncomment below

#dist_18_binding = screen_18.loc[screen_18['score'] >= 2]
#dist_18_nonbinding = screen_18.loc[screen_18['score'] <= 1]

#KDE plots for screen 15 
#binding_plot_18 = dist_18_binding['mlePeak.18'].plot.kde()
#nonbinding_plot_18 = dist_18_nonbinding['mlePeak.18'].plot.kde()

MW_stat, MW_p = sp.stats.mannwhitneyu(dist_15_binding['mlePeak.15'], dist_15_nonbinding['mlePeak.15'], alternative='two-sided')

KS_stat, KS_p = sp.stats.ks_2samp(dist_15_binding['mlePeak.15'], dist_15_nonbinding['mlePeak.15'])