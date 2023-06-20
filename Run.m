clear all

pipelinefile = 'pipe_centerline.xlsx';
wrapperfile = 'pipe_wrapper.xlsx';

cornersx = [0 1 0.7 2.3 2 3];
cornersy = [0 0 1 1 0 0];
cornersz = [0 0 1 -1 0 0];
bendradii = 0.3*[1 1 1 1];
bendpoints = 200;

pipeline = GenPipeline(pipelinefile, cornersx, cornersy, cornersz, bendradii, bendpoints);

PipeRadius = 0.1;
TurnsPerMeter = 4;
Overlap = 0.1;
Resolution = 100;
PlotAngle = pi/3;
RefVector = [0; 0; 1];

output = WrapPipe(pipeline, PipeRadius, TurnsPerMeter, Overlap, Resolution, PlotAngle, RefVector);
writetable(output, wrapperfile);