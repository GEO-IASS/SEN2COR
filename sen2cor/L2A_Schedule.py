#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

###This module starts the processes, taking the information of how many processes are alive

import os
import fnmatch
import string
from psutil import cpu_count, virtual_memory
from time import sleep
from L2A_ProcessTile import L2A_ProcessTile, SUCCESS, FAILURE
from multiprocessing import Lock, Queue, active_children

class L2A_Schedule():
    def __init__(self, config, L2A_TILES):        
        self._tiles = L2A_TILES
        self._maxProc = config.nrProcs
        self._res = config.resolution        
        self._targetDir = config.L2A_UP_DIR
        from L2A_Logger import LogQueueReader
        self._queue = Queue()
        self._msgQueue = Queue()
        log_queue_reader = LogQueueReader(self._queue)
        log_queue_reader.start()

    def sync(self):
        S2A_mask = '*S2A_*GRANULE*'
        procs = []
        for tile in self._tiles:
            if(fnmatch.fnmatch(tile, S2A_mask) == False):
                continue

            # the max processes from configuration:
            maxProc = self._maxProc
            nrProc = maxProc
        
            while True:
                # the virtual memory in percentage:
                vMemPerc = virtual_memory()[2]
                if (vMemPerc > 90.0) and (nrProc > 1):
                    nrProc -= 1
                elif (vMemPerc < 70.0) and (nrProc < maxProc):
                    nrProc += 1
                    
                if (len(procs) < nrProc):
                    p = L2A_ProcessTile(tile, self._res, self._msgQueue, self._queue)
                    try:
                        p.start()
                        procs.append(p)
                        sleep(2)
                        break
                    except:
                        continue
                else:
                    sleep(2)
                    for p in procs:
                        if (p.is_alive() == False):
                            procs.remove(p)
        
        while active_children():
            sleep(1)
            
        for t in procs:
            t.join()

        if len(procs) == 0:
            return FAILURE
        
        ret = SUCCESS
        for q in procs:
            q = self._msgQueue.get()
            if q > 0:
                ret = FAILURE
        return ret

