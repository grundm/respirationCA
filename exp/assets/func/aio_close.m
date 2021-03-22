function aio_events = aio_close(aio)
% aio_events = aio_close(aio) stops and removes analog input or output 
% objects and returns the logged events.
% 
% Author:           Martin Grund
% Last update:      October 12, 2015

%% Stop analog I/O objects

stop(aio);

%% Display events

aio_events = showdaqevents(aio);

%% Remove objects from the engine

delete(aio);