"""
A class that creates a text-based progress bar instance.

Licensed under the PSF License
source: http://code.activestate.com/recipes/473899-progress-meter/
credits:
2006-02-16 vishnubob: original code
2010-08-01 hroest: added ETA code

Here is a silly example of its usage:

import progress
import time
import random

total = 1000
p = progress.ProgressMeter(total=total)

while total > 0:
    cnt = random.randint(1, 25)
    p.update(cnt)
    total -= cnt
    time.sleep(random.random())


Here is an example of its output:

[------------------------->                                   ] 41%  821.2/sec
"""
import time, sys, math

class ProgressMeter(object):
    ESC = chr(27)
    def __init__(self, **kw):
        # What time do we start tracking our progress from?
        self.timestamp = kw.get('timestamp', time.time())
        # What kind of unit are we tracking?
        self.unit = str(kw.get('unit', ''))
        # Number of units to process
        self.total = int(kw.get('total', 100))
        # Number of units already processed
        self.count = int(kw.get('count', 0))
        # Refresh rate in seconds
        self.rate_refresh = float(kw.get('rate_refresh', .5))
        # Number of ticks in meter
        self.meter_ticks = int(kw.get('ticks', 60))
        self.meter_division = float(self.total) / self.meter_ticks
        self.meter_value = int(self.count / self.meter_division)
        self.last_update = None
        self.rate_history_idx = 0
        self.rate_history_len = 10
        self.rate_history = [None] * self.rate_history_len
        self.rate_current = 0.0
        self.last_refresh = 0
        self._cursor = False
        self.reset_cursor()

    def reset_cursor(self, first=False):
        if self._cursor:
            sys.stdout.write(self.ESC + '[u')
        self._cursor = True
        sys.stdout.write(self.ESC + '[s')

    def update(self, count, **kw):
        now = time.time()
        # Caclulate rate of progress
        rate = 0.0
        # Add count to Total
        self.count += count
        self.count = min(self.count, self.total)
        if self.last_update:
            delta = now - float(self.last_update)
            if delta:
                rate = count / delta
            else:
                rate = count
            self.rate_history[self.rate_history_idx] = rate
            self.rate_history_idx += 1
            self.rate_history_idx %= self.rate_history_len
            cnt = 0
            total = 0.0
            # Average rate history
            for rate in self.rate_history:
                if rate == None:
                    continue
                cnt += 1
                total += rate
            rate = total / cnt
        else: self.first_update = now - 0.0001
        self.rate_current = rate
        self.rate_average = self.count / (now - self.first_update)
        self.last_update = now
        # Device Total by meter division
        value = int(self.count / self.meter_division)
        if value > self.meter_value:
            self.meter_value = value
        if self.last_refresh:
            if (now - self.last_refresh) > self.rate_refresh or \
                (self.count >= self.total):
                    self.refresh()
        else:
            self.refresh()

    def nice_timestring(self, eta):
        if eta > 3600: return "%dh %dm" % ( eta // 3600, 
                               (( eta - 3600*(eta // 3600)) / 60)  )
        elif eta > 60: return "%dm %ds" % ( eta // 60, 
                               (eta - 60*(eta // 60)) )
        else: return "%ds" % eta 

    def get_meter(self, **kw):
        bar = '-' * self.meter_value
        pad = ' ' * (self.meter_ticks - self.meter_value)
        perc = (float(self.count) / self.total) * 100
        eta = (self.total - self.count ) /self.rate_average
        eta = self.nice_timestring( eta )
        return '[%s>%s] %d%%  %.1f %s/sec (eta %s)' % (bar, pad, perc, 
                   self.rate_average, self.unit,  eta )

    def refresh(self, **kw):
        # Clear line
        sys.stdout.write(self.ESC + '[2K')
        self.reset_cursor()
        sys.stdout.write(self.get_meter(**kw))
        # Are we finished?
        if self.count >= self.total:
            timediff = time.time() - self.first_update
            timediff = self.nice_timestring( timediff )
            sys.stdout.write('\n')
            sys.stdout.write('It took %s\n' % timediff )
        sys.stdout.flush()
        # Timestamp
        self.last_refresh = time.time()
