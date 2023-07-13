# PipeWrapper
Design a wrapper for an arbitrary pipe

A straight strip of material (e.g. tape) can be wrapped helically around a straight length of pipe.

This is not the case for pipes with bends.

This calculates a shape of strip that can be wrapped around a given pipe with arbitrary bends.

MATLAB code and a dumb-translated Python version with a matplotlib "GUI"

## Sample Output
### MATLAB
For the following input

```
cornersx = [0 1 0.7 2.3 2 3];
cornersy = [0 0 1 1 0 0];
cornersz = [0 0 1 -1 0 0];
bendradii = 0.3*[1 1 1 1];
bendpoints = 200;

PipeRadius = 0.1;
TurnsPerMeter = 4;
Overlap = 0.1;
Resolution = 100;
PlotAngle = pi/3;
RefVector = [0; 0; 1];
```

![3dPipe1](https://github.com/RandomVertebrate/PipeWrapper/assets/54997017/db6b5773-ef2f-4224-9859-83a27ca14921)

![2dWrapper1](https://github.com/RandomVertebrate/PipeWrapper/assets/54997017/1d43492b-01ae-42f5-bcad-bb7cc1aa9147)

### Python
![PythonOut](https://github.com/RandomVertebrate/PipeWrapper/assets/54997017/cf40502e-f6f9-4d9e-ab77-b8645b25a3d1)

### Experiment with Paper
![Experiment](https://github.com/RandomVertebrate/PipeWrapper/assets/54997017/103e2fa0-face-4ccf-988d-4d15f6255566)

