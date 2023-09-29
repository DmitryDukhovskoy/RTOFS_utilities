"""
  Relaxation utility functions
"""
def getCoord(self):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(self.img)
    cid = fig.canvas.mpl_connect('button_press_event', self.__onclick__)
    return self.point

