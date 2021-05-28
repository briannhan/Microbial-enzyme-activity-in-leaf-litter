'''Code copied and modified from following link:
https://www.geeksforgeeks.org/box-plot-in-python-using-matplotlib/
'''
# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# Creating dataset
np.random.seed(10)
data_1 = np.random.normal(100, 10, 200)
data_2 = np.random.normal(90, 20, 200)
data_3 = np.random.normal(80, 30, 200)
data_4 = np.random.normal(70, 40, 200)
data = [data_1, data_2, data_3, data_4]

fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist=True, meanline=True,
				notch=False)

colors = ['#0000FF', '#00FF00',
		'#FFFF00', '#FF00FF']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set(hatch="\\")

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
	whisker.set(color ='#8B008B',
				linewidth = 1.5,
				linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color='#8B008B', linewidth=2)

# changing color and linewidth of
# medians
for median in bp['medians']:
 	median.set(color ='red',
 			linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
	flier.set(marker ='D',
			color ='#e7298a',
			alpha = 0.5)
	
# x-axis labels
ax.set_xticklabels(['data_1', 'data_2',
					'data_3', 'data_4'])

# Adding title
plt.title("Customized box plot")
plt.legend(bp["boxes"], ["hello", "darkness", "my", "old"])

# Removing top axes and right axes
# ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
	
# show plot
# plt.show(bp)
plt.show()
