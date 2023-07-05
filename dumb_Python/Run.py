from PipeWrapSRC import *
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider, TextBox

# "GUI" for the pipe wrapper

cornersx = []
cornersy = []
cornersz = []
bendradii = []

PipeRadius = 1
TurnsPerMeter = 1
Overlap = 1
Resolution = 100
PlotAngle = 0
RefVector = [1, 1, 1]
OutputFile = "Wrapper_output.xlsx"
bendpoints = 100

current_x = 0;
current_y = 0;
current_z = 0;
current_bendrad = 0;
prev_bendrad = 0;

# Main figure (window)
fig1 = plt.figure()

# Table showing pipe centerline (input data to GenPipeline)
tableax = fig1.add_subplot(4, 2, 2)
tableax.set_title("Pipe Centerline Data")
tableax.axis("off")
tableax.table(transpose([[None],[None],[None], [None]]), colLabels = ["x", "y", "z", "Bend Radius"], loc = "center")

# Subplot for 3D view of wrapped pipe
ax3d = fig1.add_subplot(2, 2, 3, projection = "3d")
ax3d.set_title("Wrapped Pipe")

# Subplot for 2D view of unwrapped wrapper
ax2d = fig1.add_subplot(2, 2, 4)
ax2d.set_title("Wrapper")

# Update pipe centerline input data table
def updatetable():
    tableax.clear()
    tableax.axis("off")
    if len(cornersx) >= 2:
        tableax.table(transpose([cornersx, cornersy, cornersz, append(append([0], bendradii),[0])]), colLabels = ["x", "y", "z", "Bend Radius"], loc = "center")
    elif len(cornersx) == 1:
        tableax.table(transpose([cornersx, cornersy, cornersz, [0]]), colLabels = ["x", "y", "z", "bend radius"], loc = "center")
    tableax.set_title("Pipe Centerline Data")
    return

# GUI CALLBACK FUNCTIONS
# GENERATE / GEN button
def generate(event):
    if len(cornersx) < 2:
        return
    Profile = GenPipeline(cornersx, cornersy, cornersz, bendradii, bendpoints)
    WrapPipe(Profile, PipeRadius, TurnsPerMeter, Overlap, Resolution, PlotAngle, RefVector, OutputFile, ax2d, ax3d)
    return

# GENERATE / GEN button with output in new windows
def generate_newfig(event):
    if len(cornersx) < 2:
        return
    Profile = GenPipeline(cornersx, cornersy, cornersz, bendradii, bendpoints)
    WrapPipe(Profile, PipeRadius, TurnsPerMeter, Overlap, Resolution, PlotAngle, RefVector, OutputFile)
    return

# Pipe centerline x-coordinate input (update current_x)
def xinput(value):
    global current_x
    current_x = float(value)
    return

# Pipe centerline y-coordinate input (update current_y)
def yinput(value):
    global current_y
    current_y = float(value)
    return

# Pipe centerline z-coordinate input (update current_z)
def zinput(value):
    global current_z
    current_z = float(value)
    return

# Pipe centerline bend radius input (update current_bendrad)
def brinput(value):
    global current_bendrad
    current_bendrad = float(value)
    return

# Pipe radius / insulation radius input (update PipeRadius)
def prinput(value):
    global PipeRadius
    PipeRadius = float(value)

# Turns per meter input (update TurnsPerMeter)
def tpminput(value):
    global TurnsPerMeter
    TurnsPerMeter = float(value)

# Overlap input (update Overlap)
def ovinput(value):
    global Overlap
    Overlap = float(value)

# Resolution input (update Resolution)
def resinput(value):
    global Resolution
    Resolution = int(value)

# Filename input for output excel file (update OutputFile)
def fileinput(value):
    global OutputFile
    OutputFile = value

# Bend resolution input for pipe centerline generation (update bendpoints)
def bpinput(value):
    global bendpoints
    bendpoints = int(value)

# Plot Angle input for output 2D profile plotting (update PlotAngle)
def planginput(value):
    global PlotAngle
    PlotAngle = float(value)*pi/180

# RefVector x-component input (update RefVector[0])
def refxinput(value):
    global RefVector
    RefVector[0] = float(value)

# RefVector y-component input (update RefVector[1])
def refyinput(value):
    global RefVector
    RefVector[1] = float(value)

# RefVector z-component input (update RefVector[2])
def refzinput(value):
    global RefVector
    RefVector[2] = float(value)

# Add point to the pipe centerline generation input data
def addpt(event):
    global cornersx, cornersy, cornersz, bendradii
    cornersx = append(cornersx, current_x)
    cornersy = append(cornersy, current_y)
    cornersz = append(cornersz, current_z)
    if len(cornersx) > 2:
        bendradii = append(bendradii, current_bendrad)
    updatetable()

# Delete point from the pipe centerline generation input data
def delpt(event):
    global cornersx, cornersy, cornersz, bendradii
    cornersx = cornersx[:-1]
    cornersy = cornersy[:-1]
    cornersz = cornersz[:-1]
    bendradii = bendradii[:-1]
    updatetable()

# DEFINING UI ELEMENTS (INPUT TEXT BOXES, BUTTONS, SLIDERS)
# Pipe Radius input box
prinax = fig1.add_subplot(20, 6, 3)
prinbox = TextBox(prinax, "Insulation Radius:")
prinbox.on_submit(prinput)

# Turns Per Meter input box
tpminax = fig1.add_subplot(20, 6, 9)
tpminbox = TextBox(tpminax, "Turns per Unit Length:")
tpminbox.on_submit(tpminput)

# Overlap input box
ovinax = fig1.add_subplot(20, 6, 15)
ovinbox = TextBox(ovinax, "Overlap:")
ovinbox.on_submit(ovinput)

# Resolution input box
resinax = fig1.add_subplot(20, 6, 21)
resinbox = TextBox(resinax, "Resolution:", initial = 100)
resinbox.on_submit(resinput)

# Plot angle input slider
planginax = fig1.add_subplot(20, 6, 27)
planginsdr = Slider(planginax, "Plot Angle:", -90, 90, valinit = 0)
planginsdr.on_changed(planginput)

# Bend resolution input slider
bpinax = fig1.add_subplot(20, 6, 33)
bpinsdr = Slider(bpinax, "Bend Segments:", 2, 1000, valinit = 100)
bpinsdr.on_changed(bpinput)

# Output filename input box
fileinax = fig1.add_subplot(20, 6, 39)
fileinbox = TextBox(fileinax, "Output Excel File:")
fileinbox.on_submit(fileinput)

# Reference vector x-component input box
refxinax = fig1.add_subplot(20, 12, 87)
refxinbox = TextBox(refxinax, "Ref Vector:", initial = 1)
refxinbox.on_submit(refxinput)

# Reference vector y-component input box
refyinax = fig1.add_subplot(20, 12, 88)
refyinbox = TextBox(refyinax, "", initial = 1)
refyinbox.on_submit(refyinput)

# Reference vector z-component input box
refzinax = fig1.add_subplot(20, 12, 89)
refzinbox = TextBox(refzinax, "", initial = 1)
refzinbox.on_submit(refzinput)

# Pipe centerline new point x-coordinate input box
xinax = fig1.add_subplot(20, 6, 46)
xinbox = TextBox(xinax, "x:")
xinbox.on_submit(xinput)

# Pipe centerline new point y-coordinate input box
yinax = fig1.add_subplot(20, 6, 47)
yinbox = TextBox(yinax, "y:")
yinbox.on_submit(yinput)

# Pipe centerline new point z-coordinate input box
zinax = fig1.add_subplot(20, 6, 48)
zinbox = TextBox(zinax, "z:")
zinbox.on_submit(zinput)

# Bend radius input box
brinax = fig1.add_subplot(20, 6, 42)
brinbox = TextBox(brinax, "Bend Radius:")
brinbox.on_submit(brinput)

# Add point button
addbtnax = fig1.add_subplot(20, 4, 35)
addbtn = Button(addbtnax, "Add Point")
addbtn.on_clicked(addpt)

# Delete point button
delbtnax = fig1.add_subplot(20, 4, 36)
delbtn = Button(delbtnax, "Delete Point")
delbtn.on_clicked(delpt)

# GENERATE button (same window)
genbtnax = fig1.add_subplot(20, 4, 33)
genbtn = Button(genbtnax, "GEN")
genbtn.on_clicked(generate)

# GENERATE button (new window)
gennewbtnax = fig1.add_subplot(20, 4, 34)
gennewbtn = Button(gennewbtnax, "GEN (New Fig)")
gennewbtn.on_clicked(generate_newfig)

plt.subplots_adjust(hspace = 0.3)
plt.show()
