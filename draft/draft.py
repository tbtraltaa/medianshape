import numpy as np
import matplotlib.pyplot as plt

symbols = [u'\u2192']     # this helps you get the symbol 

x = np.arange(10.)
y = np.exp(-x/2.)

plt.figure()
for i, symbol in enumerate(symbols):
    y2 = np.exp(-x/2.)
    plt.plot(x, y, 'o')              # plotting with field circles
    plt.plot(x, y2, 'g')             # plotting with green line 
    for x0, y0 in zip(x, y2):
        marker = "$%s$" % symbol
        plt.plot(x0, y0, color='b', marker='>')
        plt.plot(x0, y0+10, color='b', marker='<')
plt.show()
