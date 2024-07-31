# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:52:37 2023

@author: Andres Aldana
"""
import time
import datetime


class Cronometer:

    def tick(self):
        self.ti = time.time()

    def tock(self,verbose=False):
        self.to = time.time()
        self.elapsed_seconds = self.to-self.ti
        ft = str(datetime.timedelta(seconds=self.elapsed_seconds))
        self.elapsed = ft
        if verbose:
            print("Time elapsed ",ft)


    def elapsed_milliseconds(self):
        return self.elapsed_seconds*1000


    def format_seconds(self):
        # Calculate hours, minutes, and seconds
        hours, remainder = divmod(self.elapsed_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)

        # Format and return the result
        return "{:02}:{:02}:{:02}".format(int(hours), int(minutes), int(seconds))
