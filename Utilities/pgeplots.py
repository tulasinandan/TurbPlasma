#!/usr/bin/env python
def pgeplt(rc):
    import pyqtgraph as pg
    from pyqtgraph.Qt import QtGui, QtCore
    
    rc.loadenergies()
    bd=pg.mkPen(width=2,color=(200, 200, 255), style=QtCore.Qt.DotLine)
    
    plotWidget = pg.plot(title="Change in energies for "+rc.dirname,labels={'left':'dE','bottom':'twci'})
    plotWidget.addLegend()
    plotWidget.plot(rc.t,rc.eges-rc.eges[0], pen=bd, name='dE_{tot}')
    plotWidget.plot(rc.t,rc.eb  -rc.eb[0]  , pen=20, name='dE_{b}  ')
    plotWidget.plot(rc.t,rc.eip -rc.eip[0] , pen=30, name='dE_{ip} ')
    plotWidget.plot(rc.t,rc.eep -rc.eep[0] , pen=40, name='dE_{ep} ')
    plotWidget.plot(rc.t,rc.eif -rc.eif[0] , pen=50, name='dE_{if} ')
    plotWidget.plot(rc.t,rc.eef -rc.eef[0] , pen=60, name='dE_{ef} ')
    plotWidget.plot(rc.t,rc.ee  -rc.ee[0]  , pen=70, name='dE_{ee} ') 
    QtGui.QApplication.instance().exec_()

if __name__=="__main__":
    from TurbPlasma.Utilities.subs import create_object
    rc=create_object()
    pgeplt(rc)
