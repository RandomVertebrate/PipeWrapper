clear all

pipelinefile = 'pipe_centerline.xlsx';
wrapperfile = 'pipe_wrapper.xlsx';

cornersx = [0 1 0.6 2.4 2 3];
cornersy = [0 0 1 1 0 0];
cornersz = [0 0 0 0 0 0];
bendradii = 0.2*[1 1 1 1];
bendpoints = 200;

GenPipeline(pipelinefile, cornersx, cornersy, cornersz, bendradii, bendpoints);

PipeRadius = 0.15;
TurnsPerMeter = 3;
Overlap = 0.1;
Resolution = 0.01;
PlotAngle = pi/3;
normvec = [0; 0; 1];

output = WrapPipe(pipelinefile, PipeRadius, TurnsPerMeter, Overlap, Resolution, PlotAngle, normvec)
writetable(output, wrapperfile);