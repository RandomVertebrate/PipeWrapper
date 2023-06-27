from PipeWrapSRC import *
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider, TextBox

cornersx = []
cornersy = []
cornersz = []
bendradii = []

PipeRadius = 1
TurnsPerMeter = 1
Overlap = 1
Resolution = 10
PlotAngle = 0
RefVector = [0, 0, 1]
OutputFile = "Wrapper_output.xlsx"
bendpoints = 100

current_x = 0;
current_y = 0;
current_z = 0;
current_bendrad = 0;
prev_bendrad = 0;

fig1 = plt.figure()

tableax = fig1.add_subplot(4, 2, 2)
tableax.set_title("Pipe Centerline Data")
tableax.axis("off")
tableax.table(transpose([[None],[None],[None], [None]]), colLabels = ["x", "y", "z", "Bend Radius"], loc = "center")

def generate(event):
    if len(cornersx) < 2:
        return
    Profile = GenPipeline(cornersx, cornersy, cornersz, bendradii, bendpoints)
    WrapPipe(Profile, PipeRadius, TurnsPerMeter, Overlap, Resolution, PlotAngle, RefVector, OutputFile, fig1)
    return

def updatetable():
    tableax.clear()
    tableax.axis("off")
    if len(cornersx) >= 2:
        tableax.table(transpose([cornersx, cornersy, cornersz, append(append([0], bendradii),[0])]), colLabels = ["x", "y", "z", "Bend Radius"], loc = "center")
    elif len(cornersx) == 1:
        tableax.table(transpose([cornersx, cornersy, cornersz, [0]]), colLabels = ["x", "y", "z", "bend radius"], loc = "center")
    tableax.set_title("Pipe Centerline Data")
    return

def xinput(value):
    global current_x
    current_x = float(value)
    return

def yinput(value):
    global current_y
    current_y = float(value)
    return

def zinput(value):
    global current_z
    current_z = float(value)
    return

def brinput(value):
    global current_bendrad
    current_bendrad = float(value)
    return

def prinput(value):
    global PipeRadius
    PipeRadius = float(value)

def tpminput(value):
    global TurnsPerMeter
    TurnsPerMeter = float(value)

def ovinput(value):
    global Overlap
    Overlap = float(value)

def resinput(value):
    global Resolution
    Resolution = int(value)

def fileinput(value):
    global OutputFile
    OutputFile = value

def bpinput(value):
    global bendpoints
    bendpoints = int(value)

def planginput(value):
    global PlotAngle
    PlotAngle = float(value)*pi/180

def refxinput(value):
    global RefVector
    RefVector[0] = float(value)

def refyinput(value):
    global RefVector
    RefVector[1] = float(value)

def refzinput(value):
    global RefVector
    RefVector[2] = float(value)

def addpt(event):
    global cornersx, cornersy, cornersz, bendradii
    cornersx = append(cornersx, current_x)
    cornersy = append(cornersy, current_y)
    cornersz = append(cornersz, current_z)
    if len(cornersx) > 2:
        bendradii = append(bendradii, current_bendrad)
    updatetable()

def delpt(event):
    global cornersx, cornersy, cornersz, bendradii
    cornersx = cornersx[:-1]
    cornersy = cornersy[:-1]
    cornersz = cornersz[:-1]
    bendradii = bendradii[:-1]
    updatetable()

prinax = fig1.add_subplot(20, 6, 3)
prinbox = TextBox(prinax, "Insulation Radius:")
prinbox.on_submit(prinput)

tpminax = fig1.add_subplot(20, 6, 9)
tpminbox = TextBox(tpminax, "Turns per Unit Length:")
tpminbox.on_submit(tpminput)

ovinax = fig1.add_subplot(20, 6, 15)
ovinbox = TextBox(ovinax, "Overlap:")
ovinbox.on_submit(ovinput)

resinax = fig1.add_subplot(20, 6, 21)
resinbox = TextBox(resinax, "Resolution:")
resinbox.on_submit(resinput)

planginax = fig1.add_subplot(20, 6, 27)
planginsdr = Slider(planginax, "Plot Angle:", -90, 90)
planginsdr.on_changed(planginput)

bpinax = fig1.add_subplot(20, 6, 33)
bpinsdr = Slider(bpinax, "Bend Segments:", 2, 200)
bpinsdr.on_changed(bpinput)

fileinax = fig1.add_subplot(20, 6, 39)
fileinbox = TextBox(fileinax, "Output Excel File:")
fileinbox.on_submit(fileinput)

refxinax = fig1.add_subplot(20, 12, 87)
refxinbox = TextBox(refxinax, "Ref Vector:")
refxinbox.on_submit(refxinput)

refyinax = fig1.add_subplot(20, 12, 88)
refyinbox = TextBox(refyinax, "")
refyinbox.on_submit(refyinput)

refzinax = fig1.add_subplot(20, 12, 89)
refzinbox = TextBox(refzinax, "")
refzinbox.on_submit(refzinput)

xinax = fig1.add_subplot(20, 6, 46)
xinbox = TextBox(xinax, "x:")
xinbox.on_submit(xinput)

yinax = fig1.add_subplot(20, 6, 47)
yinbox = TextBox(yinax, "y:")
yinbox.on_submit(yinput)

zinax = fig1.add_subplot(20, 6, 48)
zinbox = TextBox(zinax, "z:")
zinbox.on_submit(zinput)

brinax = fig1.add_subplot(20, 6, 42)
brinbox = TextBox(brinax, "Bend Radius:")
brinbox.on_submit(brinput)

addbtnax = fig1.add_subplot(20, 4, 35)
addbtn = Button(addbtnax, "Add Point")
addbtn.on_clicked(addpt)

delbtnax = fig1.add_subplot(20, 4, 36)
delbtn = Button(delbtnax, "Delete Point")
delbtn.on_clicked(delpt)

genbtnax = fig1.add_subplot(20, 2, 17)
genbtn = Button(genbtnax, "GENERATE")
genbtn.on_clicked(generate)

plt.subplots_adjust(hspace = 0.3)
plt.show()
