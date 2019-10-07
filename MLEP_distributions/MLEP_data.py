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
all_known_plot_15 = screen_15['mlePeak.15'].plot.kde()

#Statisical tests of binding vs. nonbinding (a)
MW_stat_a, MW_p_a = sp.stats.mannwhitneyu(dist_15_binding['mlePeak.15'], dist_15_nonbinding['mlePeak.15'], alternative='two-sided')
KS_stat_a, KS_p_a = sp.stats.ks_2samp(dist_15_binding['mlePeak.15'], dist_15_nonbinding['mlePeak.15'])

#KS_p_a = 3.168474975481327e-27
#KS_stat_a = 0.3146807109940751
#MW_p_a = 3.1892154629370733e-23
#MW_stat_a = 394413.0

#Statistical tests of binding vs. all known (b)
MW_stat_b, MW_p_b = sp.stats.mannwhitneyu(dist_15_binding['mlePeak.15'], screen_15['mlePeak.15'], alternative='two-sided')
KS_stat_b, KS_p_b = sp.stats.ks_2samp(dist_15_binding['mlePeak.15'], screen_15['mlePeak.15'])

#KS_p_b = 4.210325923173611e-18
#KS_stat_b = 0.250130821559393
#MW_p_b = 7.69044204590104e-16
#MW_stat_b = 471207.0
