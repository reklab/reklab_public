function [] = cseqPlot(C,Offset,Height,Fs,ShowMsgs, edgeColorFlag)
%PLOTPATTERNS Plots the sequence of events P in a cseq File 
%	[] = plotPatterns(P,Offset,Height,Fs,ShowMsgs) plots
%       the patterns P in the current axes.
%
%   INPUT
%   C is an N-by-1 vector with the patterns to be plotted
%       at each sample.
%   Offset is a scalar indicating the y-axis offset at
%       which the patterns will be plotted (default = 0).
%   Height is a scalar indicating the height of the
%       patterns bar (default = 1).
%   Fs is a scalar value with the sampling
%       frequency (default = 50Hz).
%   ShowMsgs is a flag indicating if messages should
%       be sent to the standard output (default = false).
%  EdgeColorFlag[false]
% - true - set edge color same
%                 - false - set edge color to default 
%
%   OUTPUTsignal
%   N/A
%
%   EXAMPLE
%   plotPatterns(P);
%
%   VERSION HISTORY
%   2015_04_29 - Created by: Carlos A. Robles-Rubio (CARR).
%
%   REFERENCES
%   [1] .
%
%
%Copyright (c) 2015-2016, Carlos Alejandro Robles Rubio, Karen A. Brown, and Robert E. Kearney, 
%McGill University
%All rights reserved.
% 
%Redistribution and use in source and binary forms, with or without modification, are 
%permitted provided that the following conditions are met:
% 
%1. Redistributions of source code must retain the above copyright notice, this list of 
%   conditions and the following disclaimer.
% 
%2. Redistributions in binary form must reproduce the above copyright notice, this list of 
%   conditions and the following disclaimer in the documentation and/or other materials 
%   provided with the distribution.
% 
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
%EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
%MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
%COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
%EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
%HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    if ~exist('Offset','var') || isempty(Offset)
        Offset=0;
    end
    if ~exist('Height','var') || isempty(Height)
        Height=1;
    end
    if ~exist('Fs','var') || isempty(Fs)
        Fs=50;
    end
    if ~exist('ShowMsgs','var') || isempty(ShowMsgs)
        ShowMsgs=false;
    end
    if ~exist('edgeColorFlag','var') || isempty(edgeColorFlag)
        edgeColorFlag=false;
    end
    E=cseq2eseq(C); 
    eseqPlot(E,Offset,Height,Fs,ShowMsgs, edgeColorFlag)
end