function lpt = dio_setup(adr1,adr2,direction)
% dio_setup creates an digital input object for the parallel port and
% controls the directionality of the data port. It returns a structure
% "lpt" that contains decimal addresses:
%   - data, control and status port (range of "adr1) 
%   - extended control register (range of "adr2")
% and functions:
%
% The parallel port must be an extended capabilities port (ECP).
% 
% Input:
% adr1          - parallel data port address range 1 in hexadecimal (e.g., '3000')
% adr2          - parallel data port address range 2 in hexadecimal (e.g., '3008')
% direction     - parallel data port mode ('uni' or 'bi')
%
% This function builds on io32.dll and inpout32.dll (for installation
% details see http://apps.usd.edu/coglab/psyc770/IO32.html)
% Now version for Matlab 32-bit on Windows 64-bit:
% http://apps.usd.edu/coglab/psyc770/IO32on64.html
%
% Administrator privileges are required.
%
% The Data Acquisition Toolbox does not support to control the 
% bidirectional mode of the parallel port and hence a modification of the 
% source code of the parallel adaptor would be necessary
% (cd([matlabroot '\toolbox\daq\daq\src\mwparallel\'])).
%
% Author:           Martin Grund
% Last update:      October 13, 2016

%% Settings
lpt.IRQ_bit = 5; % interrupt request bit (IRQ) of control port
lpt.bidir_bit = 6; % bidirectional bit of control port
% see http://retired.beyondlogic.org/spp/parallel.htm#5

lpt.mode_bits = 6:8; % of extended control register
lpt.mode_std = [0 0 0];
lpt.mode_byte = [1 0 0];

%% Addresses
%  http://retired.beyondlogic.org/ecp/ecp.htm#9
lpt.data_adr = hex2dec(adr1);
lpt.status_adr = lpt.data_adr + 1; % stauts port address
lpt.ctrl_adr = lpt.data_adr + 2; % control port address

lpt.ecr = hex2dec(adr2) + 2; % extended control register (ECR)

%% Subfunctions
lpt.get = @(dio,adr) uint8(io32(dio,adr));

lpt.set_mode = @(lpt,mode) io32(lpt.dio,lpt.ecr,set_bits(lpt.get(lpt.dio,lpt.ecr),lpt.mode_bits,mode));

lpt.set_bibit = @(lpt,bibit) io32(lpt.dio,lpt.ctrl_adr,bitset(lpt.get(lpt.dio,lpt.ctrl_adr),lpt.bidir_bit,bibit));

lpt.get_bibit = @(lpt) bitget(lpt.get(lpt.dio,lpt.ctrl_adr),lpt.bidir_bit);

lpt.set_IRQ = @(lpt,IRQbit) io32(lpt.dio,lpt.ctrl_adr,bitset(lpt.get(lpt.dio,lpt.ctrl_adr),lpt.IRQ_bit,IRQbit));

%% Create IO32 interface
clear io32;
lpt.dio = io32;

%% Install inpout32.dll driver

status = io32(lpt.dio);

if status == 0
    disp('Successful installation of inpout32a.dll'); 
else
    error('Failed installation of inpout32a.dll');
end

%% Set mode and bidirectional bit
% For modes http://retired.beyondlogic.org/ecp/ecp.htm#10)

if strcmp('uni',direction)
    % Set standard mode (bits 7:5 -> 000)
    lpt.set_mode(lpt,lpt.mode_std);
    
    % NOT NECESSARY: Turn of bidirectional bit

    if lpt.get_bibit(lpt) == 1
        error('Failed setup of parallel port unidirectional mode');
    end
    
elseif strcmp('bi',direction)
    % Set byte mode (bits 7:5 -> 001)
    lpt.set_mode(lpt,lpt.mode_byte);
    
    % Set bidirectional bit
    lpt.set_bibit(lpt,1);
    
    if lpt.get_bibit(lpt) == 0
        error('Failed setup of parallel port bidirectional mode');
    end
    
else
    error('Input argument direction has to be string with "uni" or "bi".')
end

%% Set all data pins high
io32(lpt.dio,lpt.data_adr,255);

function A = set_bits(A,bits,values)
    for i = 1:length(bits)
        A = bitset(A,bits(i),values(i));
    end