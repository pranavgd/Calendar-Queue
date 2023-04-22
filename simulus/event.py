# FILE INFO ###################################################
# Author: Jason Liu <jasonxliu2010@gmail.com>
# Created on June 14, 2019
# Last Update: Time-stamp: <2019-07-31 04:37:33 liux>
###############################################################

"""Simulation event types and event list."""

from collections.abc import MutableMapping

from .trappable import Trappable
from .trap import Trap

__all__ = ["_Event", "_DirectEvent", "_ProcessEvent", "_EventList_", \
           "infinite_time", "minus_infinite_time"]

# two extremes of simulation time
infinite_time = float('inf')
minus_infinite_time = float('-inf')

# PQDict is PRIORITY QUEUE DICTIONARY (PYTHON RECIPE)
# Created by Nezar Abdennur
# Python recipes (4591):
# http://code.activestate.com/recipes/578643-priority-queue-dictionary/
#
# An indexed priority queue implemented in pure python as a dict-like
# class. It is a stripped-down version of pqdict. A Priority Queue
# Dictionary maps dictionary keys (dkeys) to updatable priority keys
# (pkeys).
#
# The priority queue is implemented as a binary heap, which supports:
#       O(1) access to the top priority element
#       O(log n) removal of the top priority element
#       O(log n) insertion of a new element
#
# In addition, an internal dictionary or "index" maps dictionary keys
# to the position of their entry in the heap. This index is maintained
# as the heap is manipulated. As a result, a PQ-dict also supports:
#       O(1) lookup of an arbitrary element's priority key
#       O(log n) removal of an arbitrary element
#       O(log n) updating of an arbitrary element's priority key
# PQDict is modified to be used for our event list
 
class _MinEntry_(object):
    """
    Mutable entries for a Min-PQ dictionary.

    """
    def __init__(self, dkey, pkey):
        self.dkey = dkey    #dictionary key
        self.pkey = pkey    #priority key

    def __lt__(self, other):
        return self.pkey < other.pkey

class _MaxEntry_(object):
    """
    Mutable entries for a Max-PQ dictionary.

    """
    def __init__(self, dkey, pkey):
        self.dkey = dkey
        self.pkey = pkey

    def __lt__(self, other):
        return self.pkey > other.pkey

#PQDICT removed, CalQ or Sorted Dict inserted
class CalQueue():
    ''' This is a dictionary of Calender queue that maintains the position of the events inserted in
    and popped out of the Calender queue'''
    IDdict={}
    '''In Calender Queue A Hash is defined as int(time of event/width of bucket)%(Number of Buckets) '''
    def __init__(self , nBuckets ,t_width,startPrio,queuesize):
        '''m_buckets is the Array of Buckets . Buckets themselves  are basically 
        Arrays that will store  events with same Hash in a increasing order i.e, the Event with
        smallest time will be stored in index 0'''

        #print("nBuckets is ", nBuckets)
        #print("twidth is ", t_width)
        #print("startprio is ", startPrio)
        #print("queuesize ", queuesize)
        self.m_buckets=[]
        #qsize is number of events currently in calender queue
        self.m_qsize=queuesize
        for i in range (int(nBuckets)):
            self.m_buckets.append([])
        #m_nbuckets are number of buckets 
        self.m_nbuckets = int(nBuckets)
        #mwidth is width of single bucket
        self.mwidth = t_width
        #lastprio is the priority of Event which is just dequeued from the Calender Queue
        self.lastPrio=startPrio
        #Bucket top is the Priority of Bucket from which event has just been Dequeued
        self.mbucket_top=(startPrio / t_width + 1) * t_width
        index = int(startPrio/t_width)
        #mlastBucket is the bucket number from which last element is dequeued
        self.mlastBucket = index%int(nBuckets)
        
    
    '''Do Insert Function inserts the Event into calender queue based of the Hash of event 
    Within a bucket it will insert the event at position after the Event which is scheduled just before the Event
    to be inserted'''
    def DoInsert(self,entry):
        i = int(entry.time / self.mwidth)
        i = i % int(self.m_nbuckets)
        #print("i is ",i)
        for j in range(int(len(self.m_buckets[i]))):
            if self.m_buckets[i][j].time > entry.time:
                 
                self.m_buckets[i].insert(j,entry)
                CalQueue.IDdict[id(entry)]=[i,j]
                return
        CalQueue.IDdict[id(entry)]=[i,len(self.m_buckets[i])-1]
        self.m_buckets[i].append(entry)
        #print("Out of loop ")
        
        return
    def Insert(self,entry):
        self.DoInsert(entry)
        self.m_qsize=self.m_qsize+1
        #print("Current Queue Size: ",self.m_qsize,"Current Element: ",entry.time)
        self.ResizeUp(self.m_qsize)
        
        return
    '''Peek will return Pointer of Event which is Scheduled at the Earliest Time'''
    def peek(self):
        minele = float('inf')
    
        for j in range(self.m_nbuckets):
            if len(self.m_buckets[j])>0 and self.m_buckets[j][0].time < minele:
                minele = self.m_buckets[j][0].time
                minInd = j

        return self.m_buckets[minInd][0]
    '''Pop Item Will Dequeue the event Scheduled ealiest and Update the rest of bucket events'''
    def popitem(self):
        i = self.mlastBucket
        if(self.m_qsize==0) : return None
        while True:
            if (len(self.m_buckets[i])!=0 and self.m_buckets[i][0].time < self.mbucket_top):
                if id(self.m_buckets[i][0]) in CalQueue.IDdict.keys():
                    CalQueue.IDdict.pop(id(self.m_buckets[i][0]))
                temp = self.m_buckets[i].pop(0)
                self.mlastBucket = i
                self.lastPrio = temp.time
                self.m_qsize = self.m_qsize-1
                #remove the first ele out of ith bucket
                
                return temp
            else:
                i=i+1
                if i==self.m_nbuckets:
                    i=0
                self.mbucket_top=self.mbucket_top+ self.mwidth
                if i==self.mlastBucket:
                    break
        #inprio = 0
        for start in range(0,self.m_nbuckets):
            if len(self.m_buckets[start])!=0:
                self.mlastBucket = start
                self.lastPrio = self.m_buckets[start][0].time
                minprio = self.m_buckets[start][0].time
                minbucket = start
                break
        for i in range (start+1, self.m_nbuckets):
            if len(self.m_buckets[i])!=0:
                if self.m_buckets[i][0].time < minprio:
                    self.mlastBucket = i
                    self.lastPrio = self.m_buckets[i][0].time
                    minprio = self.m_buckets[i][0].time
                    minbucket = i
        if id(self.m_buckets[minbucket][0]) in CalQueue.IDdict.keys():
            CalQueue.IDdict.pop(id(self.m_buckets[minbucket][0]))
        foo = self.m_buckets[minbucket].pop(0)
        n = int(self.lastPrio / self.mwidth)
        self.mbucket_top = (n + 1.5) * self.mwidth
        self.m_qsize = self.m_qsize -1
        self.ResizeDown(self.m_qsize)
        return foo
        
    '''Deletes The Event which is Scheduled to happen Else throws error'''
    def delete_item(self, entry):
        i = (int(entry.time/self.mwidth))% self.m_nbuckets
        #bucket no. of entry
        # i = self.mbucket's First element 
        # j = self.mbucket's 
        '''
        if id(entry) in CalQueue.IDdict.keys():
            a = CalQueue.IDdict[id(entry)][0]
            b = CalQueue.IDdict[id(entry)][1]
            #print("a = ",a)
            #print("b = ",b)
            #print("Length of ath bucket ",len(self.m_buckets[a]))
            CalQueue.IDdict.pop(id(entry))
            if(b==len(self.m_buckets[a])):
                self.m_buckets[a].pop()
                self.m_qsize = self.m_qsize - 1
                self.ResizeDown(self.m_qsize)
                return
            self.m_buckets[a].pop(b)   #changed b to b-1
            self.m_qsize = self.m_qsize - 1
            self.ResizeDown(self.m_qsize)
            return
        '''
        flag = 0
        for j in range(len(self.m_buckets[i])):
            if id(entry) == id(self.m_buckets[i][j]):
                flag=1
                self.m_buckets[i].pop(j)
                CalQueue.IDdict.pop(id(entry))
                self.m_qsize = self.m_qsize - 1
                self.ResizeDown(self.m_qsize)
                break
        if flag == 0:
            raise ValueError("EventList.cancel(%s): event not found" % entry)
        else:
            return 
                
        '''
        else:
            raise ValueError("EventList.cancel(%s): event not found" % entry)'''
            
        
    
    '''According to the original document of Dr. Randy Brown Calender Queue will Change 
    number of buckets to twice when the number of Enqueued Events exceed twice the Number of Buckets
    And Halve the Number of Buckets when Number of Events Currently Scheduled is less than Half 
    the number of Buckets'''
    def ResizeUp(self,qsize):
        #print("Queue Size",qsize,"Number of buckets:",self.m_nbuckets)
        if qsize>self.m_nbuckets*2 and self.m_nbuckets < 32768:
            self.Resize(self.m_nbuckets*2,qsize)
    
    def ResizeDown(self,qsize):
        if self.m_nbuckets>=2 and qsize<self.m_nbuckets/2 :
            self.Resize(self.m_nbuckets/2,qsize)

    def Resize(self,newSize,qsize):
        newidth=self.newWidth()
        #print("newwidth",newidth)
        self.DoResize(newSize,newidth,qsize)
    ''' According to the Document of Dr.Randy Brown Width of Calender Queue is calculated by 
    taking the first 25 events Sceduled in the Calender Queue and then Calculating the Sum of 
    Time Difference of Two consequent Events and Taking twice of this sum and divinding it by number of events
    then only adding the Time Differnces if they are smaller than twice avg
    
    '''
    def newWidth(self):
        if(self.m_qsize < 2):
            return 1
        
        if(self.m_qsize <= 5):
            nsamples=self.m_qsize
        else:
            nsamples = 5 + int(self.m_qsize/10)
        if nsamples>25:
            nsamples=25
        LastBucket = self.mlastBucket
        bucketTop = self.mbucket_top
        lastprio=self.lastPrio
        samplelist=[]
        for i in range(0,nsamples):
            samplelist.append(self.popitem())
        timelist=[]
        for ele in samplelist:
            timelist.append(ele.time)
        #print("Sample List: ",timelist)
        for obj in samplelist:
            self.DoInsert(obj)
        self.mlastBucket = LastBucket
        self.mbucket_top = bucketTop
        self.lastPrio = lastprio
        
        cur = 0
        nxt = 1
        totalsep = 0
        while nxt!=len(samplelist):
            totalsep = totalsep + samplelist[nxt].time - samplelist[cur].time
            cur = cur + 1
            nxt= nxt + 1
        
        twiceavg = totalsep /(nsamples - 1)*2
        totalsep = 0
        cur = 0
        nxt = 1
        while nxt!=len(samplelist):
            diff = samplelist[nxt].time - samplelist[cur].time
            if diff<=twiceavg:
                totalsep = totalsep + diff
            cur=cur+1
            nxt=nxt+1
        totalsep = totalsep*3
        totalsep = max(totalsep, 1)
        #print("tot sep ", totalsep)
        return int(totalsep)

    def DoResize(self, newsize, newwidth,qsize):
        oldbuckets  = self.m_buckets
        oldN = self.m_nbuckets
        oldq=qsize
        #print("oldq",oldq)
        oldbuckets  = self.m_buckets
        oldN = self.m_nbuckets
        
        #print(storeself)
        #python internally stores object and its assoc attributes
        # as a dictionary as key-val pairs 
        newObj = CalQueue(newsize, newwidth, self.lastPrio,oldq)
        self.__dict__.update(newObj.__dict__)
        for i in range(oldN):
            for j in range(len(oldbuckets[i])):
                #print("doresize ",oldbuckets[i][j])
                self.DoInsert(oldbuckets[i][j])
        
        
        #setattr(self, old)
    def update_item(self, entry):
        #print("update item called ")
        for i in range(self.m_nbuckets):
            for j in range(len(self.m_buckets[i])):
                if id(entry)== id(self.m_buckets[i][j]):
                    self.m_buckets[i].pop(j)
                    CalQueue.IDdict.pop(id(entry))
                    self.Insert(entry)
                    return

    def printInfo(self):
        print("No. of buckets = ", self.m_nbuckets)
        print("No. of Event = ", self.m_qsize) 
        print("Bucket width = ", self.mwidth)
        print("Last priority = ", self.lastPrio)
        #print(" = %d", self.m_nbuckets)
        #print("No. of buckets = %d", self.m_nbuckets)
        for i in range(self.m_nbuckets):
            print("Day ", i, end=" :  ")
            for j in range(len(self.m_buckets[i])):
                print( self.m_buckets[i][j].time, end="  ")
            print()
### PQDict or Sorted Dict ends
class _Event(Trappable):
    """The base class for all simulation events."""

    def __init__(self, sim, time, name=None):
        super().__init__(sim)
        self.time = time
        self.name = name
        self.trap = None

    def __str__(self):
        return "%g: evt=%s" % \
            (self.time, self.name if self.name else id(self))

    def __lt__(self, other):
        return self.time < other.time

    def _try_wait(self):
        # if the current event is on the event list, pass onto the
        # trap created for this event; otherwise, we consider the
        # trappable already triggered
        if self._sim._eventlist.current_event(self):
            if self.trap is None:
                self.trap = Trap(self._sim)
            return self.trap._try_wait()
        else:
            return False
        
    def _cancel_wait(self):
        assert self.trap is not None
        self.trap._cancel_wait()

    def _true_trappable(self):
        # it's possible the event is no longer in the event list, in
        # which case the trap is None (meaning it's sprung)
        if self.trap is not None:
            return self.trap
        else:
            return self
        
class _DirectEvent(_Event):
    """The event type for direct event scheduling."""

    #def __init__(self, sim, time, func, params, name, repeat_intv):
    def __init__(self, sim, time, func, name, repeat_intv, usr_args, usr_kwargs):
        super().__init__(sim, time, name)
        self.func = func
        #self.params = params
        self.repeat_intv = repeat_intv
        self.args = usr_args
        self.kwargs = usr_kwargs

    def __str__(self):
        return "%g: dir_evt=%s %s" % \
            (self.time, self.name if self.name else self.func.__name__+'()',
             "(repeat=%g)"%self.repeat_intv if self.repeat_intv else "")

    def renew(self, time):
        self.time = time
        self.trap = None # trap cannot be reused
        return self

class _ProcessEvent(_Event):
    """The event type for process scheduling."""

    def __init__(self, sim, time, proc, name):
        #print("process event called")
        super().__init__(sim, time, name)
        self.proc = proc

    def __str__(self):
        return "%g: prc_evt=%s" % \
            (self.time, self.name if self.name else self.proc.func.__name__+'()')

# Modified Event List
class _EventList_(object):
    
    """An event list sorts events in timestamp order.
    An event list is a priority queue that stores and sorts simulation
    events based on the time of the events. It supports three basic
    operations: to insert a (future) event, to peek and to retrieve
    the event with the minimal timestamp.
    """
    
    def __init__(self):
        #self.pqueue = []
        self.cqueue = CalQueue(2,1,0,0)
        self.EvtDict=CalQueue.IDdict
        self.last =float('-inf')

    def __len__(self):
        return self.cqueue.m_qsize
        
    def insert(self, evt):
        #print('insert called')
        if self.last <= evt.time:
            #heapq.heappush(self.pqueue, (evt.time, id(evt), evt))
            self.cqueue.Insert(evt)
        else:
            raise ValueError("EventList.insert(%s): past event (last=%g)" %
                             (evt, self.last))

    def get_min(self):
        #print("peek called ")
        if self.cqueue.m_qsize > 0:
            #return self.pqueue[0][0] # just return the time
            e = self.cqueue.peek()
            return e.time  # just return the time
        else:
            raise IndexError("EventList.get_min() from empty list")

    def delete_min(self):
        #print("pop item called")
        if self.cqueue.m_qsize > 0:
            #assert self.last <= self.pqueue[0][0]
            #self.last = self.pqueue[0][0]
            #return heapq.heappop(self.pqueue)[-1]
            e = self.cqueue.popitem()
            #print("self", self.last, "event time ", e.time)
            assert self.last <= e.time
            self.last = e.time
            return e
        else:
            raise IndexError("EventList.delete_min() from empty list")

    def cancel(self, evt):
        #print("delete called")
        self.cqueue.delete_item(evt)

    def update(self, evt):
        #print("update called")
        if id(evt) not in self.EvtDict:
            raise ValueError("EventList.update(%s): event not found" % evt)
        if self.last <= evt.time:
            self.cqueue.update_item(evt)
            #self.cqueue.Insert(evt)
        else:
            raise ValueError("EventList.update(%s): past event (last=%g)" %
                             (evt, self.last))

    def current_event(self, evt):
        # check whether the event is current
        return id(evt) in self.EvtDict